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


#include "BoysFuncGPU.hpp"
#include "EriExchangeHessianSXPX.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_II_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x, PA_y;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            // i-i Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        (-2.0) * a_i * (
                            +QC_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S1 * a_i * a_i * (
                            +QC_0*delta[g0][g1]
                        )

                        + 4.0 * a_i * a_i * (
                            +PA_x*PA_y*QC_0
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                            +QC_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_i * (
                            -PQ[c0]*delta[g0][g1]

                            +PA_x*delta[c0][g1] + PA_y*delta[c0][g0]
                        )

                        + (-4.0) * S1 * inv_S4 * a_i * a_i * (
                            +PA_x*PA_y*PQ[c0]
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PQ[c0]*delta[g0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_i * (
                            +PA_x*PQ[g1]*QC_0

                            +PA_y*PQ[g0]*QC_0
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PA_x*PQ[c0]*PQ[g1]

                            +PA_y*PQ[c0]*PQ[g0]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F3_t[3] * (

                        (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ii_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_KK_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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




        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        (-2.0) * a_k * (
                            +QC_0*delta[g0][g1]

                            +QC_x*delta[c0][g1] + QC_y*delta[c0][g0]
                        )

                        + 2.0 * inv_S2 * a_k * a_k * (
                            +QC_0*delta[g0][g1]

                            +QC_x*delta[c0][g1] + QC_y*delta[c0][g0]
                        )

                        + 4.0 * a_k * a_k * (
                            +QC_0*QC_x*QC_y
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_k * (
                            +delta[g0][g1]*(PQ[c0] + QC_0)

                            +delta[c0][g0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[g0] + QC_x)
                        )

                        + (-4.0) * S1 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*QC_x*QC_y

                            +QC_0*(PQ[g0]*QC_y + PQ[g1]*QC_x)
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                    )

                    + F3_t[2] * (

                        2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*PQ[g0]*QC_y

                            +PQ[g1]*(PQ[c0]*QC_x + PQ[g0]*QC_0)
                        )

                    )

                    + F3_t[3] * (

                        (-4.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kk_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k], hess_kk_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_IK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
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

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        (-2.0) * a_i * (
                            +PA_x*delta[c0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_k * (
                            +PA_x*delta[c0][g1]
                        )

                        + 4.0 * a_i * a_k * (
                            +PA_x*QC_0*QC_y
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_k * (
                            +PA_x*delta[c0][g1]
                        )

                        + (-4.0) * S1 * inv_S4 * a_i * a_k * (
                            +PA_x*PQ[c0]*QC_y

                            +PA_x*PQ[g1]*QC_0
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_k * (
                            +PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_k * (
                            +PQ[g0]*QC_0*QC_y
                        )

                    )

                    + F3_t[2] * (

                        (-2.0) * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*PQ[g0]*QC_y

                            +PQ[g0]*PQ[g1]*QC_0
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PQ[c0]*PQ[g1]
                        )

                    )

                    + F3_t[3] * (

                        4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ik_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ik_xy += ERIs[y][x];
            }
        }

        // Note factor of 2 due to IK<->JL symmetry for ground state Hessian

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_IJ_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x, PB_y;
        uint32_t j_prim, j_cgto;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            // i-j Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        4.0 * a_i * a_j * (
                            +PA_x*PB_y*QC_0
                        )

                        + 2.0 * inv_S1 * a_i * a_j * (
                            +QC_0*delta[g0][g1]
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                            +QC_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_j * (
                            +PB_y*delta[c0][g0]

                            -PQ[c0]*delta[g0][g1]

                            +PA_x*delta[c0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_j * (
                            -PA_x*PB_y*PQ[c0]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_j * (
                            +QC_0*(PA_x*PQ[g1] + PB_y*PQ[g0])
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PA_x*PQ[c0]*PQ[g1]

                            +PB_y*PQ[c0]*PQ[g0]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F3_t[3] * (

                        (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ij_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_j,
            hess_ij_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_KL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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




        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_ps[displ_k + l];

            // Q_kl == Q_K_ps[displ_k + l]
            if (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        2.0 * inv_S2 * a_k * a_l * (
                            +QD_y*delta[c0][g0]

                            +QC_0*delta[g0][g1] + QC_x*delta[c0][g1]
                        )

                        + 4.0 * a_k * a_l * (
                            +QC_0*QC_x*QD_y
                        )

                        + (-2.0) * a_l * (
                            +QD_y*delta[c0][g0]
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_l * (
                            +delta[g0][g1]*(PQ[c0] + QC_0)

                            +delta[c0][g0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[g0] + QC_x)
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_l * (
                            -PQ[g0]*QC_0*QD_y - QC_x*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_l * (
                            +PQ[g1]*delta[c0][g0]
                        )

                    )

                    + F3_t[2] * (

                        4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*QC_0
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                    )

                    + F3_t[3] * (

                        (-4.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kl_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kl_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_l,
            hess_kl_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_IL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_ps[displ_k + l];

            // Q_kl == Q_K_ps[displ_k + l]
            if (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        4.0 * a_i * a_l * (
                            +PA_x*QC_0*QD_y
                        )

                        + 2.0 * inv_S2 * a_i * a_l * (
                            +PA_x*delta[c0][g1]
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_l * (
                            +PA_x*delta[c0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_l * (
                            +QD_y*delta[c0][g0]

                            +PQ[g0]*delta[c0][g1]

                            +QC_0*delta[g0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_l * (
                            -PA_x*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_l * (
                            +PQ[g0]*QC_0*QD_y
                        )

                    )

                    + F3_t[2] * (

                        (-2.0) * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PQ[g0]*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PQ[c0]*PQ[g1]
                        )

                    )

                    + F3_t[3] * (

                        4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_il_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_il_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_l,
            hess_il_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + atom_idx_l * natoms + prim_cart_ao_to_atom_inds[i],
            hess_il_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPS_JK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ss_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_x;
        uint32_t j_prim, j_cgto;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

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

            PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ss_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, rho * d2 * r2_PQ, 3, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F3_t[1] *= d2;
                F3_t[2] *= d2 * d2;
                F3_t[3] *= d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // j-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F3_t[0] * (

                        2.0 * inv_S2 * a_j * a_k * (
                            +PB_x*delta[c0][g1]
                        )

                        + 4.0 * a_j * a_k * (
                            +PB_x*QC_0*QC_y
                        )

                        + (-2.0) * a_j * (
                            +PB_x*delta[c0][g1]
                        )

                    )

                    + F3_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_j * a_k * (
                            +PB_x*delta[c0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_j * a_k * (
                            -PB_x*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + (-2.0) * S2 * inv_S4 * a_j * (
                            +PQ[g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S4 * a_j * a_k * (
                            +PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]
                        )

                        + 4.0 * S2 * inv_S4 * a_j * a_k * (
                            +PQ[g0]*QC_0*QC_y
                        )

                    )

                    + F3_t[2] * (

                        (-2.0) * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*delta[g0][g1]

                            +PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_x*PQ[c0]*PQ[g1]
                        )

                        + (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*PQ[g0]*QC_y

                            +PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F3_t[3] * (

                        4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_jk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_jk_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + atom_idx_j * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_jk_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_j,
            hess_jk_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_II_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x, PA_y;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            // i-i Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S1 * a_i * a_i * (
                            +QC_0*QD_0*delta[g0][g1]
                        )

                        + (-1.0) * inv_S2 * a_i * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + 4.0 * a_i * a_i * (
                            +PA_x*PA_y*QC_0*QD_0
                        )

                        + (-2.0) * a_i * (
                            +QC_0*QD_0*delta[g0][g1]
                        )

                        + inv_S1 * inv_S2 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_i * (
                            +PA_x*PA_y*delta[c0][d0]
                        )

                    )

                    + F4_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_i * (
                            +PA_x*PA_y*delta[c0][d0]
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                            +QC_0*QD_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_i * (
                            +PA_x*QD_0*delta[c0][g1]

                            +PA_y*QD_0*delta[c0][g0] + QC_0*(PA_x*delta[d0][g1] + PA_y*delta[d0][g0])

                            -delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +delta[c0][d0]*(PA_x*PQ[g1] + PA_y*PQ[g0])
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_i * (
                            -PA_x*PA_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_i * (
                            +QC_0*QD_0*(PA_x*PQ[g1] + PA_y*PQ[g0])
                        )

                        + S1 * inv_S2 * inv_S4 * a_i * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PA_x*PQ[c0]*delta[d0][g1]

                            -PA_y*PQ[c0]*delta[d0][g0] - PQ[d0]*(PA_x*delta[c0][g1] + PA_y*delta[c0][g0])

                            -delta[c0][d0]*(PA_x*PQ[g1] + PA_y*PQ[g0])

                            +PQ[c0]*PQ[d0]*delta[g0][g1]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            +PQ[g0]*PQ[g1]*delta[c0][d0]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            -(PA_x*PQ[g1] + PA_y*PQ[g0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[g0]*PQ[g1]*QC_0*QD_0
                        )

                        + inv_S4 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PA_x*PA_y*PQ[c0]*PQ[d0]
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PA_x*PQ[c0]*PQ[d0]*PQ[g1]

                            +PA_y*PQ[c0]*PQ[d0]*PQ[g0]
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ii_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_KK_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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




        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_k * a_k * (
                            +QC_0*QD_0*delta[g0][g1]

                            +QC_x*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0] + QD_0*delta[c0][g1]) + QC_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + 4.0 * a_k * a_k * (
                            +QC_0*QC_x*QC_y*QD_0
                        )

                        + (-2.0) * a_k * (
                            +QC_0*QD_0*delta[g0][g1]

                            +QD_0*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                        + inv_S2 * inv_S2 * a_k * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S2 * inv_S4 * a_k * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_k * (
                            +delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QC_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0))
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_k * (
                            -QC_0*QD_0*(PQ[g0]*QC_y + PQ[g1]*QC_x) - QC_x*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[d0]*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0]) + QD_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +delta[g0][g1]*(PQ[c0]*(PQ[d0] + QD_0) + PQ[d0]*QC_0)

                            +PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QC_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*QC_x*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[g0]*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0) + PQ[g1]*QC_0*(PQ[d0]*QC_x + PQ[g0]*QD_0)
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[d0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                        + (-4.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*QC_y

                            +PQ[g1]*(PQ[c0]*(PQ[d0]*QC_x + PQ[g0]*QD_0) + PQ[d0]*PQ[g0]*QC_0)
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kk_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k], hess_kk_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_IK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
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

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_i * a_k * (
                            +PA_x*QD_0*delta[c0][g1]

                            +PA_x*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + 4.0 * a_i * a_k * (
                            +PA_x*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * a_i * (
                            +PA_x*QD_0*delta[c0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_k * (
                            +PA_x*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PA_x*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 2.0 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*delta[d0][g0] + QD_0*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0])

                            +PQ[g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[c0][g1]*delta[d0][g0]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_k * (
                            -PA_x*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_k * (
                            +PQ[g0]*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[g0]*QD_0*delta[c0][g1]
                        )

                        + inv_S2 * inv_S4 * a_i * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*PQ[d0]*delta[c0][g1]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            +PQ[g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PQ[g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PQ[c0]*delta[d0][g1]

                            +PA_x*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[d0]*PQ[g0]*delta[c0][g1]
                        )

                    )

                    + F4_t[3] * (

                        (-4.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ik_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ik_xy += ERIs[y][x];
            }
        }

        // Note factor of 2 due to IK<->JL symmetry for ground state Hessian

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_IJ_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x, PB_y;
        uint32_t j_prim, j_cgto;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            // i-j Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S1 * a_i * a_j * (
                            +QC_0*QD_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_j * (
                            +PA_x*PB_y*delta[c0][d0]
                        )

                        + 4.0 * a_i * a_j * (
                            +PA_x*PB_y*QC_0*QD_0
                        )

                        + inv_S1 * inv_S2 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_j * (
                            +PA_x*PB_y*delta[c0][d0]
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                            +QC_0*QD_0*delta[g0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_j * (
                            +PA_x*QD_0*delta[c0][g1]

                            +delta[c0][d0]*(PA_x*PQ[g1] + PB_y*PQ[g0])

                            +PB_y*QD_0*delta[c0][g0] + QC_0*(PA_x*delta[d0][g1] + PB_y*delta[d0][g0])

                            -delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_j * (
                            -PA_x*PB_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_j * (
                            +QC_0*QD_0*(PA_x*PQ[g1] + PB_y*PQ[g0])
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PA_x*PQ[c0]*delta[d0][g1]

                            -PB_y*PQ[c0]*delta[d0][g0] - PQ[d0]*(PA_x*delta[c0][g1] + PB_y*delta[c0][g0])

                            -delta[c0][d0]*(PA_x*PQ[g1] + PB_y*PQ[g0])

                            +PQ[c0]*PQ[d0]*delta[g0][g1]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            +PQ[g0]*PQ[g1]*delta[c0][d0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PA_x*PB_y*PQ[c0]*PQ[d0]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            -(PA_x*PQ[g1] + PB_y*PQ[g0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[g0]*PQ[g1]*QC_0*QD_0
                        )

                        + inv_S4 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*(PA_x*PQ[g1] + PB_y*PQ[g0])
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ij_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_j,
            hess_ij_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_KL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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




        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_pp[displ_k + l];

            // Q_kl == Q_K_pp[displ_k + l]
            if (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_k * a_l * (
                            +QC_0*QD_0*delta[g0][g1]

                            +QC_x*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1] + QD_y*delta[c0][d0]) + QD_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                        + (-1.0) * inv_S2 * a_l * (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                        + 4.0 * a_k * a_l * (
                            +QC_0*QC_x*QD_0*QD_y
                        )

                        + (-2.0) * a_k * (
                            +QC_0*QC_x*delta[d0][g1]
                        )

                        + (-2.0) * a_l * (
                            +QD_0*QD_y*delta[c0][g0]
                        )

                        + inv_S2 * inv_S2 * a_k * a_l * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S2 * inv_S4 * a_k * a_l * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_l * (
                            +delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QD_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0))
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_l * (
                            -QC_0*QD_y*(PQ[d0]*QC_x + PQ[g0]*QD_0) - QC_x*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_l * (
                            +delta[c0][g0]*(PQ[d0]*QD_y + PQ[g1]*QD_0)
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                        + S1 * inv_S2 * inv_S4 * a_l * (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +delta[d0][g1]*(PQ[c0]*QC_x + PQ[g0]*QC_0)
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +delta[g0][g1]*(PQ[c0]*(PQ[d0] + QD_0) + PQ[d0]*QC_0)

                            +PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QD_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*QC_x*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*QC_0*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +PQ[c0]*PQ[g0]*delta[d0][g1]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_l * (
                            +PQ[d0]*PQ[g1]*delta[c0][g0]
                        )

                        + S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PQ[c0]*PQ[d0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) - PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kl_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kl_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_l,
            hess_kl_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_IL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        uint32_t j_prim, j_cgto;

        if (j < count_i)
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

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_pp[displ_k + l];

            // Q_kl == Q_K_pp[displ_k + l]
            if (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_i * a_l * (
                            +PA_x*QD_0*delta[c0][g1]

                            +PA_x*(QC_0*delta[d0][g1] + QD_y*delta[c0][d0])
                        )

                        + 4.0 * a_i * a_l * (
                            +PA_x*QC_0*QD_0*QD_y
                        )

                        + (-2.0) * a_i * (
                            +PA_x*QC_0*delta[d0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_l * (
                            +PA_x*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PA_x*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 2.0 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_y*delta[d0][g0] + QD_0*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QD_y*delta[c0][g0])

                            +PQ[g0]*(QC_0*delta[d0][g1] + QD_y*delta[c0][d0])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[c0][g0]*delta[d0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_l * (
                            -PA_x*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_l * (
                            +PQ[g0]*QC_0*QD_0*QD_y
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[g0]*QC_0*delta[d0][g1]
                        )

                        + inv_S2 * inv_S4 * a_i * a_l * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*PQ[c0]*delta[d0][g1]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*(QD_0*delta[g0][g1] + QD_y*delta[d0][g0]) + PQ[d0]*(QC_0*delta[g0][g1] + QD_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            +PQ[g0]*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PQ[g0]*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PQ[c0]*delta[d0][g1]

                            +PA_x*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*PQ[g0]*delta[d0][g1]
                        )

                    )

                    + F4_t[3] * (

                        (-4.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[g0]*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_il_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_il_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_l,
            hess_il_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + atom_idx_l * natoms + prim_cart_ao_to_atom_inds[i],
            hess_il_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSSPP_JK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_ss,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_ss,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_ss,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_ss[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_ss[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_x;
        uint32_t j_prim, j_cgto;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

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

            PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);



        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // j-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_j * a_k * (
                            +PB_x*QC_0*delta[d0][g1]

                            +PB_x*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + 4.0 * a_j * a_k * (
                            +PB_x*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * a_j * (
                            +PB_x*QD_0*delta[c0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_j * a_k * (
                            +PB_x*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PB_x*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 2.0 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*delta[d0][g0] + QD_0*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0])

                            +PQ[g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + (-1.0) * inv_S4 * a_j * (
                            +delta[c0][g1]*delta[d0][g0]
                        )

                        + 4.0 * S1 * inv_S4 * a_j * a_k * (
                            -PB_x*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_j * (
                            +PB_x*PQ[d0]*delta[c0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_j * a_k * (
                            +PQ[g0]*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * S2 * inv_S4 * a_j * (
                            +PQ[g0]*QD_0*delta[c0][g1]
                        )

                        + inv_S2 * inv_S4 * a_j * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S1 * inv_S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_x*PQ[c0]*delta[d0][g1]

                            +PB_x*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            +PQ[g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_x*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PQ[g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * (
                            +PQ[d0]*PQ[g0]*delta[c0][g1]
                        )

                    )

                    + F4_t[3] * (

                        4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PB_x*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_jk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_jk_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + atom_idx_j * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_jk_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_j,
            hess_jk_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_II_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x, PA_y;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            // i-i Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S1 * a_i * a_i * (
                            +PB_0*QC_0*delta[g0][g1]

                            +QC_0*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0])
                        )

                        + 4.0 * a_i * a_i * (
                            +PA_x*PA_y*PB_0*QC_0
                        )

                        + (-2.0) * a_i * (
                            +PB_0*QC_0*delta[g0][g1]
                        )

                    )

                    + F4_t[1] * (

                        2.0 * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                            +QC_0*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +QC_0*(delta[b0][g0]*(-PA_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * inv_S4 * a_i * a_i * (
                            +PA_x*PB_0*delta[c0][g1]

                            +PA_y*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -PQ[c0]*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0] + PB_0*delta[g0][g1])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[b0][c0]*delta[g0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_i * (
                            -PA_x*PA_y*PB_0*PQ[c0]
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PB_0*PQ[c0]*delta[g0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_i * (
                            +QC_0*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[b0]*QC_0*delta[g0][g1]
                        )

                        + inv_S1 * inv_S4 * a_i * a_i * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + (-2.0) * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*QC_0*delta[g0][g1]

                            +QC_0*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*delta[g0][g1]*(PB_0 - PQ[b0])

                            +PA_x*(PQ[b0]*delta[c0][g1] + PQ[g1]*delta[b0][c0]) + PA_y*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]) + PB_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])

                            +PQ[c0]*(delta[b0][g0]*(PA_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PQ[c0]*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +QC_0*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]
                        )

                    )

                    + F4_t[3] * (

                        4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PQ[c0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))
                        )

                        + 2.0 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ii_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_KK_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];


            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_k * a_k * (
                            +PB_0*QC_0*delta[g0][g1]

                            +PB_0*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                        + 4.0 * a_k * a_k * (
                            +PB_0*QC_0*QC_x*QC_y
                        )

                        + (-2.0) * a_k * (
                            +PB_0*QC_0*delta[g0][g1]

                            +PB_0*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_k * (
                            +PB_0*delta[g0][g1]*(PQ[c0] + QC_0)

                            +PB_0*(delta[c0][g0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[g0] + QC_x))
                        )

                        + (-1.0) * inv_S4 * a_k * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_k * (
                            -PB_0*(PQ[g0]*QC_0*QC_y + QC_x*(PQ[c0]*QC_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +PB_0*PQ[c0]*delta[g0][g1]

                            +PB_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_k * (
                            +PQ[b0]*QC_0*delta[g0][g1]

                            +PQ[b0]*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                        + inv_S2 * inv_S4 * a_k * a_k * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * inv_S4 * a_k * a_k * (
                            +QC_0*(PQ[b0]*delta[g0][g1] + QC_x*delta[b0][g1] + QC_y*delta[b0][g0]) + QC_x*QC_y*delta[b0][c0]

                            +PQ[b0]*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                        + 4.0 * S2 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*QC_0*QC_x*QC_y
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*PQ[c0]*delta[g0][g1]

                            +PB_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0] + QC_0) + PQ[c0]*(QC_x*delta[b0][g1] + QC_y*delta[b0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][g0] + QC_x*delta[b0][c0])

                            +PQ[b0]*(delta[c0][g0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[g0] + QC_x))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*(PQ[c0]*(PQ[g0]*QC_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*QC_0)
                        )

                        + (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*QC_x*QC_y

                            +PQ[b0]*QC_0*(PQ[g0]*QC_y + PQ[g1]*QC_x)
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[b0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                    )

                    + F4_t[3] * (

                        4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            -PB_0*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*QC_y

                            +PQ[b0]*PQ[g1]*(PQ[c0]*QC_x + PQ[g0]*QC_0)
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kk_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k], hess_kk_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_IK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        (-1.0) * inv_S1 * a_i * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_k * (
                            +PA_x*PB_0*delta[c0][g1]
                        )

                        + 4.0 * a_i * a_k * (
                            +PA_x*PB_0*QC_0*QC_y
                        )

                        + (-2.0) * a_i * (
                            +PA_x*PB_0*delta[c0][g1]
                        )

                        + inv_S1 * inv_S2 * a_i * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S1 * a_i * a_k * (
                            +QC_0*QC_y*delta[b0][g0]
                        )

                    )

                    + F4_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_k * (
                            +PA_x*PB_0*delta[c0][g1]
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*delta[b0][g0]
                        )

                        + 2.0 * inv_S4 * a_i * a_k * (
                            +delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            +PA_x*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PB_0*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0])

                            -delta[b0][g0]*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_k * (
                            -PA_x*PB_0*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + S2 * inv_S1 * inv_S4 * a_i * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PB_0*PQ[c0]*delta[g0][g1]

                            -delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            -PA_x*(PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) - PB_0*PQ[g1]*delta[c0][g0]

                            +PQ[c0]*PQ[g1]*delta[b0][g0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PB_0*PQ[c0]*PQ[g1]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + (-2.0) * S2 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[g0]*delta[c0][g1]
                        )

                        + inv_S4 * inv_S4 * a_i * a_k * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0])

                            +delta[b0][g0]*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[g0]*QC_0*QC_y
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[c0]*PQ[g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*QC_y

                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ik_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ik_xy += ERIs[y][x];
            }
        }

        // Note factor of 2 due to IK<->JL symmetry for ground state Hessian

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_IJ_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x, PB_y;
        uint32_t j_prim, j_cgto, b0;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            // i-j Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S1 * a_i * a_j * (
                            +PB_0*QC_0*delta[g0][g1]

                            +QC_0*(PA_x*delta[b0][g1] + PB_y*delta[b0][g0])
                        )

                        + 4.0 * a_i * a_j * (
                            +PA_x*PB_0*PB_y*QC_0
                        )

                        + (-2.0) * a_i * (
                            +PA_x*QC_0*delta[b0][g1]
                        )

                    )

                    + F4_t[1] * (

                        2.0 * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                            +QC_0*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +QC_0*(delta[b0][g0]*(-PB_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * inv_S4 * a_i * a_j * (
                            +PA_x*PB_0*delta[c0][g1]

                            +PB_y*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -PQ[c0]*(PA_x*delta[b0][g1] + PB_0*delta[g0][g1] + PB_y*delta[b0][g0])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[b0][g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_j * (
                            -PA_x*PB_0*PB_y*PQ[c0]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_j * (
                            +QC_0*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[g0]*QC_0*delta[b0][g1]
                        )

                        + inv_S1 * inv_S4 * a_i * a_j * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*PQ[c0]*delta[b0][g1]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + (-2.0) * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*QC_0*delta[g0][g1]

                            +QC_0*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*delta[g0][g1]*(PB_0 - PQ[b0])

                            +PA_x*(PQ[b0]*delta[c0][g1] + PQ[g1]*delta[b0][c0]) + PB_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]) + PB_y*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0])

                            +PQ[c0]*(delta[b0][g0]*(PB_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PQ[c0]*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +QC_0*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*PQ[g0]*delta[b0][g1]
                        )

                    )

                    + F4_t[3] * (

                        4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PQ[c0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))
                        )

                        + 2.0 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ij_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_j,
            hess_ij_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_KL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];


            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_ps[displ_k + l];

            // Q_kl == Q_K_ps[displ_k + l]
            if (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S2 * a_k * a_l * (
                            +PB_0*QC_0*delta[g0][g1]

                            +PB_0*(QC_x*delta[c0][g1] + QD_y*delta[c0][g0])
                        )

                        + 4.0 * a_k * a_l * (
                            +PB_0*QC_0*QC_x*QD_y
                        )

                        + (-2.0) * a_l * (
                            +PB_0*QD_y*delta[c0][g0]
                        )

                    )

                    + F4_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_l * (
                            +PB_0*delta[g0][g1]*(PQ[c0] + QC_0)

                            +PB_0*(delta[c0][g0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[g0] + QC_x))
                        )

                        + 2.0 * inv_S4 * a_k * a_l * (
                            +QC_0*QC_x*delta[b0][g1] + QD_y*(PQ[b0]*delta[c0][g0] + QC_0*delta[b0][g0] + QC_x*delta[b0][c0])

                            +PQ[b0]*(QC_0*delta[g0][g1] + QC_x*delta[c0][g1])
                        )

                        + (-1.0) * inv_S4 * a_l * (
                            +delta[b0][g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_l * (
                            -PB_0*(PQ[g0]*QC_0*QD_y + QC_x*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_l * (
                            +PB_0*PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S2 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*QC_0*QC_x*QD_y
                        )

                        + (-2.0) * S2 * inv_S4 * a_l * (
                            +PQ[b0]*QD_y*delta[c0][g0]
                        )

                        + inv_S2 * inv_S4 * a_k * a_l * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                    )

                    + F4_t[2] * (

                        (-1.0) * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*PQ[c0]*delta[g0][g1]

                            +PB_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0] + QC_0) + PQ[c0]*(QC_x*delta[b0][g1] + QD_y*delta[b0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QD_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][g0] + QC_x*delta[b0][c0])

                            +PQ[b0]*(delta[c0][g0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[g0] + QC_x))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*(PQ[c0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PQ[b0]*(PQ[g0]*QC_0*QD_y + QC_x*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_l * (
                            +PQ[b0]*PQ[g1]*delta[c0][g0]
                        )

                    )

                    + F4_t[3] * (

                        4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PB_0*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*(PQ[c0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                    )

                    + F4_t[4] * (

                        (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kl_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kl_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_l,
            hess_kl_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_IL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_ps[displ_k + l];

            // Q_kl == Q_K_ps[displ_k + l]
            if (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        2.0 * inv_S1 * a_i * a_l * (
                            +QC_0*QD_y*delta[b0][g0]
                        )

                        + 2.0 * inv_S2 * a_i * a_l * (
                            +PA_x*PB_0*delta[c0][g1]
                        )

                        + 4.0 * a_i * a_l * (
                            +PA_x*PB_0*QC_0*QD_y
                        )

                        + inv_S1 * inv_S2 * a_i * a_l * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_l * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_l * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_l * (
                            +PA_x*PB_0*delta[c0][g1]
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_y*delta[b0][g0]
                        )

                        + 2.0 * inv_S4 * a_i * a_l * (
                            +delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            +PA_x*(QC_0*delta[b0][g1] + QD_y*delta[b0][c0]) + PB_0*(QC_0*delta[g0][g1] + QD_y*delta[c0][g0])

                            -delta[b0][g0]*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_l * (
                            -PA_x*PB_0*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_y*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PB_0*PQ[c0]*delta[g0][g1]

                            -delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            -PA_x*(PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) - PB_0*PQ[g1]*delta[c0][g0]

                            +PQ[c0]*PQ[g1]*delta[b0][g0]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*QD_y*delta[c0][g0]

                            +delta[b0][g0]*(PQ[c0]*QD_y + PQ[g1]*QC_0)

                            +PQ[g0]*QD_y*delta[b0][c0] + QC_0*(PQ[b0]*delta[g0][g1] + PQ[g0]*delta[b0][g1])

                            +PQ[b0]*PQ[g0]*delta[c0][g1]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PB_0*PQ[c0]*PQ[g1]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[g0]*QC_0*QD_y
                        )

                        + inv_S4 * inv_S4 * a_i * a_l * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[c0]*PQ[g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PQ[b0]*PQ[g0]*(PQ[c0]*QD_y + PQ[g1]*QC_0)
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_il_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_il_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_l,
            hess_il_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + atom_idx_l * natoms + prim_cart_ao_to_atom_inds[i],
            hess_il_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPS_JK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    ps_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_ps,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_ps,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_ps,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_ps,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_ps,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_ps[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_ps[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PB_x;
        uint32_t j_prim, j_cgto, b0;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_ps[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_ps[displ_k + l] * ps_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_ps[displ_k + l];

            const auto l_prim = D_inds_K_ps[displ_k + l];

            const auto l_cgto = s_prim_aoinds[l_prim];

            const auto a_l = s_prim_info[l_prim + s_prim_count * 0];

            const double r_l[3] = {s_prim_info[l_prim + s_prim_count * 2],
                                   s_prim_info[l_prim + s_prim_count * 3],
                                   s_prim_info[l_prim + s_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_ps[displ_k + l];


            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            double d2 = 1.0;

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

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // j-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F4_t[0] * (

                        (-1.0) * inv_S1 * a_j * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S2 * a_j * a_k * (
                            +PB_0*PB_x*delta[c0][g1]
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 4.0 * a_j * a_k * (
                            +PB_0*PB_x*QC_0*QC_y
                        )

                        + (-2.0) * a_j * (
                            +PB_0*PB_x*delta[c0][g1]
                        )

                        + (-2.0) * a_k * (
                            +QC_0*QC_y*delta[b0][g0]
                        )

                        + inv_S1 * inv_S2 * a_j * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S1 * a_j * a_k * (
                            +QC_0*QC_y*delta[b0][g0]
                        )

                        + (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                    )

                    + F4_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_j * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_j * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_j * a_k * (
                            +PB_0*PB_x*delta[c0][g1]
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*delta[b0][g0]
                        )

                        + 2.0 * inv_S4 * a_j * a_k * (
                            +delta[c0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])

                            +PB_0*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PB_x*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0])

                            -delta[b0][g0]*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * inv_S4 * a_j * a_k * (
                            -PB_0*PB_x*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_j * (
                            +delta[c0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + S2 * inv_S1 * inv_S4 * a_j * (
                            +delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +delta[b0][g0]*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                    )

                    + F4_t[2] * (

                        2.0 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PB_0*PQ[c0]*delta[g0][g1]

                            -delta[c0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])

                            -PB_x*PQ[c0]*delta[b0][g1] - PQ[g1]*(PB_0*delta[c0][g0] + PB_x*delta[b0][c0])

                            +PQ[c0]*PQ[g1]*delta[b0][g0]
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_0*PB_x*PQ[c0]*PQ[g1]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +PQ[c0]*PQ[g1]*delta[b0][g0]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            -(PB_0*PQ[g0] + PB_x*PQ[b0])*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + (-2.0) * S2 * S2 * inv_S4 * inv_S4 * a_j * (
                            +PQ[b0]*PQ[g0]*delta[c0][g1]
                        )

                        + inv_S4 * inv_S4 * a_j * a_k * (
                            +delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0])

                            +delta[b0][g0]*(PQ[c0]*QC_y + PQ[g1]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[g0]*QC_0*QC_y
                        )

                    )

                    + F4_t[3] * (

                        (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[c0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*delta[c0][g1] + PQ[c0]*delta[b0][g1] + PQ[g1]*delta[b0][c0]) + PQ[g1]*(PQ[b0]*delta[c0][g0] + PQ[c0]*delta[b0][g0])
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[c0]*PQ[g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*QC_y

                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0
                        )

                    )

                    + F4_t[4] * (

                        4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_jk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_jk_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + atom_idx_j * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_jk_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_j,
            hess_jk_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_II_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x, PA_y;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            // i-i Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S1 * inv_S2 * a_i * a_i * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0])
                        )

                        + 2.0 * inv_S1 * a_i * a_i * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]

                            +QC_0*QD_0*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0])
                        )

                        + 2.0 * inv_S2 * a_i * a_i * (
                            +PA_x*PA_y*PB_0*delta[c0][d0]
                        )

                        + (-1.0) * inv_S2 * a_i * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]
                        )

                        + 4.0 * a_i * a_i * (
                            +PA_x*PA_y*PB_0*QC_0*QD_0
                        )

                        + (-2.0) * a_i * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]
                        )

                    )

                    + F5_t[1] * (

                        inv_S1 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +delta[c0][d0]*(delta[b0][g0]*(-PA_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_i * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_i * (
                            +PA_x*PA_y*PB_0*delta[c0][d0]
                        )

                        + S1 * inv_S2 * inv_S4 * a_i * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]
                        )

                        + 2.0 * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                            +QC_0*QD_0*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +QC_0*QD_0*(delta[b0][g0]*(-PA_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])

                            +PA_x*(PA_y*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]) + PB_0*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1])) + PA_y*PB_0*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            -(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0] + PB_0*delta[g0][g1])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[g0][g1]*(PQ[b0]*delta[c0][d0] + QC_0*delta[b0][d0] + QD_0*delta[b0][c0])
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_i * (
                            -PA_x*PA_y*PB_0*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PB_0*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_i * (
                            +QC_0*QD_0*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[b0]*QC_0*QD_0*delta[g0][g1]
                        )

                    )

                    + F5_t[2] * (

                        (-1.0) * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                        + inv_S4 * inv_S4 * a_i * a_i * (
                            +delta[c0][d0]*delta[g0][g1]*(PB_0 - PQ[b0])

                            +PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1]) + PA_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0]) + PB_0*(delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) - PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])

                            +delta[c0][d0]*(delta[b0][g0]*(PA_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))
                        )

                        + 2.0 * S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PA_x*PB_0*PQ[c0]*delta[d0][g1]

                            -PA_y*(PA_x*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0]) + PB_0*PQ[c0]*delta[d0][g0]) - PB_0*PQ[d0]*(PA_x*delta[c0][g1] + PA_y*delta[c0][g0])

                            -delta[c0][d0]*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])

                            +PQ[c0]*PQ[d0]*(PA_x*delta[b0][g1] + PA_y*delta[b0][g0] + PB_0*delta[g0][g1])
                        )

                        + (-2.0) * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*QC_0*QD_0*delta[g0][g1]

                            +QC_0*QD_0*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +delta[g0][g1]*(PB_0 - PQ[b0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +delta[c0][d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))

                            +PA_x*(PQ[b0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) + PA_y*(PQ[b0]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) + PB_0*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))

                            +(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(delta[b0][g0]*(PA_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PA_x*PA_y*PB_0*PQ[c0]*PQ[d0]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_i * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            -(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                            +QC_0*QD_0*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))
                        )

                        + S1 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*delta[c0][d0]*delta[g0][g1]

                            +delta[g0][g1]*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0])
                        )

                    )

                    + F5_t[3] * (

                        2.0 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]*(-PB_0 + PQ[b0])

                            -PA_x*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1]) + PQ[g1]*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0])) - PA_y*(PQ[b0]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0]) + PQ[g0]*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0])) - PB_0*(PQ[c0]*(PQ[g0]*delta[d0][g1] + PQ[g1]*delta[d0][g0]) + PQ[d0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]))

                            -delta[c0][d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))

                            +PQ[c0]*PQ[d0]*(delta[b0][g0]*(-PA_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[b0]*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PQ[g0]*PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])

                            +PQ[b0]*PQ[g0]*PQ[g1]*delta[c0][d0]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[d0]*(PA_x*(PA_y*PQ[b0] + PB_0*PQ[g1]) + PA_y*PB_0*PQ[g0])
                        )

                        + (-2.0) * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            -(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))
                        )

                        + 4.0 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0*QD_0
                        )

                        + S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                    )

                    + F5_t[4] * (

                        (-2.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                        + 4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[c0]*PQ[d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PA_y*PQ[g0]))
                        )

                        + 4.0 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            -PQ[b0]*PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                    )

                    + F5_t[5] * (

                        4.0 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ii_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_KK_0(double*         hess_xy,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];


            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S2 * inv_S2 * a_k * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + 2.0 * inv_S2 * a_k * a_k * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]

                            +PB_0*(QC_x*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0] + QD_0*delta[c0][g1]) + QC_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + 4.0 * a_k * a_k * (
                            +PB_0*QC_0*QC_x*QC_y*QD_0
                        )

                        + (-2.0) * a_k * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]

                            +PB_0*QD_0*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                    )

                    + F5_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S2 * inv_S4 * a_k * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + inv_S2 * inv_S4 * a_k * a_k * (
                            +QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])

                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QC_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + QC_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_k * (
                            +PB_0*delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PB_0*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QC_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0)))
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + 2.0 * inv_S4 * a_k * a_k * (
                            +QC_0*QD_0*(PQ[b0]*delta[g0][g1] + QC_x*delta[b0][g1] + QC_y*delta[b0][g0]) + QC_x*QC_y*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +PQ[b0]*(QC_x*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0] + QD_0*delta[c0][g1]) + QC_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))
                        )

                        + (-1.0) * inv_S4 * a_k * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]) + QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])

                            +delta[b0][d0]*(QC_0*delta[g0][g1] + QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_k * (
                            -PB_0*(QC_0*QD_0*(PQ[g0]*QC_y + PQ[g1]*QC_x) + QC_x*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +PB_0*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PB_0*(PQ[d0]*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0]) + QD_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]))
                        )

                        + 4.0 * S2 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*QC_0*QC_x*QC_y*QD_0
                        )

                        + (-2.0) * S2 * inv_S4 * a_k * (
                            +PQ[b0]*QC_0*QD_0*delta[g0][g1]

                            +PQ[b0]*QD_0*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0])
                        )

                    )

                    + F5_t[2] * (

                        S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            -2*PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +(-PQ[c0] - QC_0)*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + (-PQ[d0] - QD_0)*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + (-PQ[g0] - QC_x)*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + (-PQ[g1] - QC_y)*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*delta[g0][g1]*(PQ[c0]*(PQ[d0] + QD_0) + PQ[d0]*QC_0)

                            +PB_0*(PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QC_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0]))
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PQ[b0]*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QC_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0)))

                            +QC_0*(PQ[d0]*(QC_x*delta[b0][g1] + QC_y*delta[b0][g0]) + PQ[g0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[g1]*(QC_x*delta[b0][d0] + QD_0*delta[b0][g0])) + QC_x*(PQ[c0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + delta[b0][c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0)) + QC_y*QD_0*(PQ[c0]*delta[b0][g0] + PQ[g0]*delta[b0][c0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*(PQ[c0]*QC_x*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[g0]*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0) + PQ[g1]*QC_0*(PQ[d0]*QC_x + PQ[g0]*QD_0))
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PB_0*PQ[d0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * a_k * (
                            -PQ[b0]*(QC_0*QD_0*(PQ[g0]*QC_y + PQ[g1]*QC_x) + QC_x*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[b0]*(PQ[d0]*(QC_x*delta[c0][g1] + QC_y*delta[c0][g0]) + QD_0*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0]))
                        )

                        + S1 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +delta[b0][d0]*(PQ[c0]*delta[g0][g1] + PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])

                            +PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                    )

                    + F5_t[3] * (

                        (-2.0) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PB_0*(PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0]))
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[c0]*(PQ[b0]*delta[g0][g1]*(PQ[d0] + QD_0) + PQ[d0]*(QC_x*delta[b0][g1] + QC_y*delta[b0][g0]) + PQ[g0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[g1]*(QC_x*delta[b0][d0] + QD_0*delta[b0][g0])) + PQ[d0]*(QC_0*(PQ[b0]*delta[g0][g1] + PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0]) + delta[b0][c0]*(PQ[g0]*QC_y + PQ[g1]*QC_x)) + PQ[g0]*PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +PQ[b0]*(PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QC_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0]))
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            -PB_0*(PQ[c0]*PQ[d0]*(PQ[g0]*QC_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*(PQ[c0]*QC_x*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[g0]*QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0) + PQ[g1]*QC_0*(PQ[d0]*QC_x + PQ[g0]*QD_0))
                        )

                        + (-2.0) * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[b0]*PQ[d0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])
                        )

                        + S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                    )

                    + F5_t[4] * (

                        (-2.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                        + 4.0 * S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PB_0*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                        + (-4.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*QC_y

                            +PQ[b0]*PQ[g1]*(PQ[c0]*(PQ[d0]*QC_x + PQ[g0]*QD_0) + PQ[d0]*PQ[g0]*QC_0)
                        )

                    )

                    + F5_t[5] * (

                        4.0 * S1 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kk_xy += ERIs[y][x];
            }
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k], hess_kk_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_IK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S1 * inv_S2 * a_i * a_k * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]

                            +delta[b0][g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + 2.0 * inv_S1 * a_i * a_k * (
                            +QC_0*QC_y*QD_0*delta[b0][g0]
                        )

                        + (-1.0) * inv_S1 * a_i * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_k * (
                            +PA_x*PB_0*QC_0*delta[d0][g1]

                            +PA_x*PB_0*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + 4.0 * a_i * a_k * (
                            +PA_x*PB_0*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * a_i * (
                            +PA_x*PB_0*QD_0*delta[c0][g1]
                        )

                    )

                    + F5_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_k * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]

                            +delta[b0][g0]*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + inv_S2 * inv_S4 * a_i * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))

                            +PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_k * (
                            +PA_x*PB_0*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PA_x*PB_0*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*QD_0*delta[b0][g0]
                        )

                        + S2 * inv_S1 * inv_S4 * a_i * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S4 * a_i * a_k * (
                            +QC_0*delta[d0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])

                            +QC_0*(PA_x*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PB_0*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1])) + QC_y*QD_0*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -delta[b0][g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + inv_S4 * a_i * (
                            -PB_0*delta[c0][g1]*delta[d0][g0]

                            -PA_x*delta[b0][d0]*delta[c0][g1]

                            +PQ[d0]*delta[b0][g0]*delta[c0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_k * (
                            -PA_x*PB_0*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*PB_0*PQ[d0]*delta[c0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*QD_0*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +QD_0*delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                    )

                    + F5_t[2] * (

                        S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            +delta[b0][g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + inv_S4 * inv_S4 * a_i * a_k * (
                            +QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g1]*delta[c0][g0])

                            +delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))

                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g1]*delta[d0][g0]) + QC_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0])
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PB_0*PQ[c0]*delta[d0][g1]

                            +PA_x*PB_0*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            -delta[d0][g1]*(PQ[c0] + QC_0)*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            -PA_x*(PQ[c0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[d0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) - PB_0*(PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))

                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +QC_0*QC_y*(PQ[b0]*delta[d0][g0] + PQ[g0]*delta[b0][d0]) + QD_0*(PQ[b0]*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]))

                            +delta[b0][g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))

                            +PQ[b0]*PQ[g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + (-1.0) * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*delta[c0][g1]*delta[d0][g0]

                            +delta[c0][g1]*(PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PA_x*PB_0*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[d0]*delta[c0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[g0]*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * S2 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[g0]*QD_0*delta[c0][g1]
                        )

                    )

                    + F5_t[3] * (

                        (-1.0) * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])

                            +PQ[c0]*(PA_x*(PQ[d0]*delta[b0][g1] + PQ[g1]*delta[b0][d0]) + PB_0*PQ[g1]*delta[d0][g0]) + PQ[d0]*PQ[g1]*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -PQ[c0]*PQ[d0]*PQ[g1]*delta[b0][g0]
                        )

                        + (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*(PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PQ[g0]*(PQ[c0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[d0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]))

                            +PQ[b0]*PQ[g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PA_x*PB_0*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PQ[b0]*PQ[g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[d0]*PQ[g0]*delta[c0][g1]
                        )

                    )

                    + F5_t[4] * (

                        4.0 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            -PQ[c0]*PQ[d0]*PQ[g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + 4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                    )

                    + F5_t[5] * (

                        (-4.0) * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ik_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ik_xy += ERIs[y][x];
            }
        }

        // Note factor of 2 due to IK<->JL symmetry for ground state Hessian

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ik_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_IJ_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x, PB_y;
        uint32_t j_prim, j_cgto, b0;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
            PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            // i-j Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S1 * inv_S2 * a_i * a_j * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PA_x*delta[b0][g1] + PB_y*delta[b0][g0])
                        )

                        + 2.0 * inv_S1 * a_i * a_j * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]

                            +QC_0*QD_0*(PA_x*delta[b0][g1] + PB_y*delta[b0][g0])
                        )

                        + 2.0 * inv_S2 * a_i * a_j * (
                            +PA_x*PB_0*PB_y*delta[c0][d0]
                        )

                        + (-1.0) * inv_S2 * a_i * (
                            +PA_x*delta[b0][g1]*delta[c0][d0]
                        )

                        + 4.0 * a_i * a_j * (
                            +PA_x*PB_0*PB_y*QC_0*QD_0
                        )

                        + (-2.0) * a_i * (
                            +PA_x*QC_0*QD_0*delta[b0][g1]
                        )

                    )

                    + F5_t[1] * (

                        inv_S1 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +delta[c0][d0]*(delta[b0][g0]*(-PB_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                        + (-1.0) * inv_S2 * inv_S4 * a_i * a_j * (
                            +PB_0*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PA_x*delta[b0][g1] + PB_y*delta[b0][g0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_j * (
                            +PA_x*PB_0*PB_y*delta[c0][d0]
                        )

                        + 2.0 * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                            +QC_0*QD_0*delta[g0][g1]*(-PB_0 + PQ[b0])

                            +QC_0*QD_0*(delta[b0][g0]*(-PB_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])

                            +PA_x*(PB_0*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PB_y*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) + PB_0*PB_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                            -(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PA_x*delta[b0][g1] + PB_0*delta[g0][g1] + PB_y*delta[b0][g0])
                        )

                        + (-1.0) * inv_S4 * a_i * (
                            +delta[b0][g1]*(PQ[g0]*delta[c0][d0] + QC_0*delta[d0][g0] + QD_0*delta[c0][g0])
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_j * (
                            -PA_x*PB_0*PB_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*delta[b0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_j * (
                            +QC_0*QD_0*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +PQ[g0]*QC_0*QD_0*delta[b0][g1]
                        )

                        + S1 * inv_S2 * inv_S4 * a_i * (
                            +PA_x*delta[b0][g1]*delta[c0][d0]
                        )

                    )

                    + F5_t[2] * (

                        (-1.0) * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*delta[c0][d0]*delta[g0][g1]

                            +delta[c0][d0]*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                        + inv_S4 * inv_S4 * a_i * a_j * (
                            +delta[c0][d0]*delta[g0][g1]*(PB_0 - PQ[b0])

                            +PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1]) + PB_0*(delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]) + PB_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0])

                            +delta[c0][d0]*(delta[b0][g0]*(PB_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))

                            -PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) - PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])
                        )

                        + 2.0 * S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PA_x*PB_0*PQ[c0]*delta[d0][g1]

                            -PB_y*PQ[c0]*(PA_x*delta[b0][d0] + PB_0*delta[d0][g0]) - PQ[d0]*(PA_x*(PB_0*delta[c0][g1] + PB_y*delta[b0][c0]) + PB_0*PB_y*delta[c0][g0])

                            -delta[c0][d0]*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])

                            +PQ[c0]*PQ[d0]*(PA_x*delta[b0][g1] + PB_0*delta[g0][g1] + PB_y*delta[b0][g0])
                        )

                        + (-2.0) * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*QC_0*QD_0*delta[g0][g1]

                            +QC_0*QD_0*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +delta[g0][g1]*(PB_0 - PQ[b0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +delta[c0][d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))

                            +PA_x*(PQ[b0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) + PB_0*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PB_y*(PQ[b0]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]))

                            +(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(delta[b0][g0]*(PB_y - PQ[g1]) + delta[b0][g1]*(PA_x - PQ[g0]))
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PA_x*PB_0*PB_y*PQ[c0]*PQ[d0]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_i * (
                            +PA_x*PQ[c0]*PQ[d0]*delta[b0][g1]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            -(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[g0]*delta[b0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                            +QC_0*QD_0*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))
                        )

                        + S1 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*delta[b0][g1]*delta[d0][g0]

                            +delta[b0][g1]*(PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0])
                        )

                    )

                    + F5_t[3] * (

                        2.0 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*delta[g0][g1]*(-PB_0 + PQ[b0])

                            -PA_x*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1]) + PQ[g1]*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0])) - PB_0*(PQ[c0]*(PQ[g0]*delta[d0][g1] + PQ[g1]*delta[d0][g0]) + PQ[d0]*(PQ[g0]*delta[c0][g1] + PQ[g1]*delta[c0][g0])) - PB_y*(PQ[b0]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0]) + PQ[g0]*(PQ[c0]*delta[b0][d0] + PQ[d0]*delta[b0][c0]))

                            -delta[c0][d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))

                            +PQ[c0]*PQ[d0]*(delta[b0][g0]*(-PB_y + PQ[g1]) + delta[b0][g1]*(-PA_x + PQ[g0]))
                        )

                        + 2.0 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                            +PQ[b0]*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PQ[g0]*PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0])

                            +PQ[b0]*PQ[g0]*PQ[g1]*delta[c0][d0]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*(PA_x*(PB_0*PQ[g1] + PB_y*PQ[b0]) + PB_0*PB_y*PQ[g0])
                        )

                        + (-2.0) * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*PQ[d0]*PQ[g0]*delta[b0][g1]
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            -(PQ[c0]*QD_0 + PQ[d0]*QC_0)*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))
                        )

                        + 4.0 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[g0]*PQ[g1]*QC_0*QD_0
                        )

                        + S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                    )

                    + F5_t[4] * (

                        (-2.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                        + 4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[c0]*PQ[d0]*(PB_0*PQ[g0]*PQ[g1] + PQ[b0]*(PA_x*PQ[g1] + PB_y*PQ[g0]))
                        )

                        + 4.0 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            -PQ[b0]*PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                        )

                    )

                    + F5_t[5] * (

                        4.0 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_ij_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_j,
            hess_ij_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_KL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];


            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_pp[displ_k + l];

            // Q_kl == Q_K_pp[displ_k + l]
            if (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_x = (a_l * inv_S2) * (r_l[g0] - r_k[g0]);
            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // k-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S2 * inv_S2 * a_k * a_l * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + 2.0 * inv_S2 * a_k * a_l * (
                            +PB_0*QC_0*QD_0*delta[g0][g1]

                            +PB_0*(QC_x*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1] + QD_y*delta[c0][d0]) + QD_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +PB_0*delta[c0][g0]*delta[d0][g1]
                        )

                        + (-1.0) * inv_S2 * a_l * (
                            +PB_0*delta[c0][g0]*delta[d0][g1]
                        )

                        + 4.0 * a_k * a_l * (
                            +PB_0*QC_0*QC_x*QD_0*QD_y
                        )

                        + (-2.0) * a_k * (
                            +PB_0*QC_0*QC_x*delta[d0][g1]
                        )

                        + (-2.0) * a_l * (
                            +PB_0*QD_0*QD_y*delta[c0][g0]
                        )

                        + (
                            +PB_0*delta[c0][g0]*delta[d0][g1]
                        )

                    )

                    + F5_t[1] * (

                        (-2.0) * S1 * inv_S2 * inv_S2 * inv_S4 * a_k * a_l * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + inv_S2 * inv_S4 * a_k * a_l * (
                            +QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + QC_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + QD_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])

                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_k * a_l * (
                            +PB_0*delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PB_0*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QD_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0)))
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +PB_0*delta[c0][g0]*delta[d0][g1]
                        )

                        + S1 * inv_S2 * inv_S4 * a_l * (
                            +PB_0*delta[c0][g0]*delta[d0][g1]
                        )

                        + 2.0 * inv_S4 * a_k * a_l * (
                            +QC_0*QD_0*(PQ[b0]*delta[g0][g1] + QC_x*delta[b0][g1] + QD_y*delta[b0][g0]) + QC_x*QD_y*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +PQ[b0]*(QC_x*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1] + QD_y*delta[c0][d0]) + QD_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))
                        )

                        + (-1.0) * inv_S4 * a_k * (
                            +delta[d0][g1]*(PQ[b0]*delta[c0][g0] + QC_0*delta[b0][g0] + QC_x*delta[b0][c0])
                        )

                        + (-1.0) * inv_S4 * a_l * (
                            +delta[c0][g0]*(PQ[b0]*delta[d0][g1] + QD_0*delta[b0][g1] + QD_y*delta[b0][d0])
                        )

                        + 4.0 * S1 * inv_S4 * a_k * a_l * (
                            -PB_0*(QC_0*QD_y*(PQ[d0]*QC_x + PQ[g0]*QD_0) + QC_x*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +PB_0*delta[d0][g1]*(PQ[c0]*QC_x + PQ[g0]*QC_0)
                        )

                        + 2.0 * S1 * inv_S4 * a_l * (
                            +PB_0*delta[c0][g0]*(PQ[d0]*QD_y + PQ[g1]*QD_0)
                        )

                        + 4.0 * S2 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*QC_0*QC_x*QD_0*QD_y
                        )

                        + (-2.0) * S2 * inv_S4 * a_k * (
                            +PQ[b0]*QC_0*QC_x*delta[d0][g1]
                        )

                        + (-2.0) * S2 * inv_S4 * a_l * (
                            +PQ[b0]*QD_0*QD_y*delta[c0][g0]
                        )

                        + S2 * inv_S4 * (
                            +PQ[b0]*delta[c0][g0]*delta[d0][g1]
                        )

                    )

                    + F5_t[2] * (

                        S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])
                        )

                        + S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            -2*PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +(-PQ[c0] - QC_0)*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + (-PQ[d0] - QD_0)*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + (-PQ[g0] - QC_x)*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + (-PQ[g1] - QD_y)*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*delta[g0][g1]*(PQ[c0]*(PQ[d0] + QD_0) + PQ[d0]*QC_0)

                            +PB_0*(PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QD_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0]))
                        )

                        + (-2.0) * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*delta[g0][g1]*(PQ[d0]*QC_0 + QD_0*(PQ[c0] + QC_0))

                            +PQ[b0]*(PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + QC_x*(delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0)) + QD_y*(delta[c0][d0]*(PQ[g0] + QC_x) + delta[c0][g0]*(PQ[d0] + QD_0) + delta[d0][g0]*(PQ[c0] + QC_0)))

                            +QC_0*(PQ[d0]*(QC_x*delta[b0][g1] + QD_y*delta[b0][g0]) + PQ[g0]*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + PQ[g1]*(QC_x*delta[b0][d0] + QD_0*delta[b0][g0])) + QC_x*(PQ[c0]*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + delta[b0][c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0)) + QD_0*QD_y*(PQ[c0]*delta[b0][g0] + PQ[g0]*delta[b0][c0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*(PQ[c0]*QC_x*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*QC_0*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +PB_0*PQ[c0]*PQ[g0]*delta[d0][g1]
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_l * (
                            +PB_0*PQ[d0]*PQ[g1]*delta[c0][g0]
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PQ[b0]*(QC_0*QD_y*(PQ[d0]*QC_x + PQ[g0]*QD_0) + QC_x*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_l * (
                            +PQ[b0]*delta[c0][g0]*(PQ[d0]*QD_y + PQ[g1]*QD_0)
                        )

                        + S1 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*delta[c0][g0]*delta[d0][g1]

                            +delta[d0][g1]*(PQ[c0]*delta[b0][g0] + PQ[g0]*delta[b0][c0])
                        )

                        + S1 * inv_S4 * inv_S4 * a_l * (
                            +PQ[b0]*delta[c0][g0]*delta[d0][g1]

                            +delta[c0][g0]*(PQ[d0]*delta[b0][g1] + PQ[g1]*delta[b0][d0])
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*delta[d0][g1]*(PQ[c0]*QC_x + PQ[g0]*QC_0)
                        )

                    )

                    + F5_t[3] * (

                        (-2.0) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PB_0*(PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0]))
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[c0]*(PQ[b0]*delta[g0][g1]*(PQ[d0] + QD_0) + PQ[d0]*(QC_x*delta[b0][g1] + QD_y*delta[b0][g0]) + PQ[g0]*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + PQ[g1]*(QC_x*delta[b0][d0] + QD_0*delta[b0][g0])) + PQ[d0]*(QC_0*(PQ[b0]*delta[g0][g1] + PQ[g0]*delta[b0][g1] + PQ[g1]*delta[b0][g0]) + delta[b0][c0]*(PQ[g0]*QD_y + PQ[g1]*QC_x)) + PQ[g0]*PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])

                            +PQ[b0]*(PQ[g1]*(QC_0*delta[d0][g0] + QC_x*delta[c0][d0] + QD_0*delta[c0][g0]) + delta[c0][g1]*(PQ[d0]*(PQ[g0] + QC_x) + PQ[g0]*QD_0) + delta[d0][g1]*(PQ[c0]*(PQ[g0] + QC_x) + PQ[g0]*QC_0) + (PQ[g1] + QD_y)*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0] + PQ[g0]*delta[c0][d0]))
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PB_0*(PQ[c0]*PQ[d0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*(PQ[c0]*QC_x*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*QC_0*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + (-2.0) * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*delta[d0][g1]
                        )

                        + (-2.0) * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_l * (
                            +PQ[b0]*PQ[d0]*PQ[g1]*delta[c0][g0]
                        )

                        + S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                    )

                    + F5_t[4] * (

                        (-2.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                        + 4.0 * S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PB_0*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            -PQ[b0]*(PQ[c0]*PQ[d0]*(PQ[g0]*QD_y + PQ[g1]*QC_x) + PQ[g0]*PQ[g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                    )

                    + F5_t[5] * (

                        4.0 * S1 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_kl_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_kl_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_l,
            hess_kl_xy * ik_factor_D * 2.0 * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_IL_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_l,
                                const uint32_t* atom_ids_l_displs,
                                const uint32_t* atom_ids_l_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_l = 0; atom_idx_l < natoms; atom_idx_l++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_l_displ = atom_ids_l_displs[k * natoms + atom_idx_l];
    const auto atom_l_count = atom_ids_l_counts[k * natoms + atom_idx_l];

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PA_x;
        uint32_t j_prim, j_cgto, b0;

        if (j < count_i)
        {
            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (atom_l_count + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const auto l_raw = n * TILE_DIM_X_K + threadIdx.x;

            if ((j >= count_i) || (l_raw >= atom_l_count))
            {
                break;
            }

            // note non-coalesced read wrt l
            const auto l = atom_ids_l[displ_k + atom_l_displ + l_raw];

            // const auto Q_kl = Q_K_pp[displ_k + l];

            // Q_kl == Q_K_pp[displ_k + l]
            if (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold)
            {
                break;
            }

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

            // i-l Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S1 * inv_S2 * a_i * a_l * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]

                            +delta[b0][g0]*(QC_0*delta[d0][g1] + QD_y*delta[c0][d0])
                        )

                        + 2.0 * inv_S1 * a_i * a_l * (
                            +QC_0*QD_0*QD_y*delta[b0][g0]
                        )

                        + (-1.0) * inv_S1 * a_i * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]
                        )

                        + 2.0 * inv_S2 * a_i * a_l * (
                            +PA_x*PB_0*QC_0*delta[d0][g1]

                            +PA_x*PB_0*(QD_0*delta[c0][g1] + QD_y*delta[c0][d0])
                        )

                        + 4.0 * a_i * a_l * (
                            +PA_x*PB_0*QC_0*QD_0*QD_y
                        )

                        + (-2.0) * a_i * (
                            +PA_x*PB_0*QC_0*delta[d0][g1]
                        )

                    )

                    + F5_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_i * a_l * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]

                            +delta[b0][g0]*(QD_0*delta[c0][g1] + QD_y*delta[c0][d0])
                        )

                        + inv_S2 * inv_S4 * a_i * a_l * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))

                            +PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_i * a_l * (
                            +PA_x*PB_0*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PA_x*PB_0*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_0*QD_y*delta[b0][g0]
                        )

                        + 2.0 * inv_S4 * a_i * a_l * (
                            +QC_0*delta[d0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(QD_0*delta[c0][g1] + QD_y*delta[c0][d0])

                            +QC_0*(PA_x*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + PB_0*(QD_0*delta[g0][g1] + QD_y*delta[d0][g0])) + QD_0*QD_y*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -delta[b0][g0]*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + inv_S4 * a_i * (
                            -PB_0*delta[c0][g0]*delta[d0][g1]

                            -PA_x*delta[b0][c0]*delta[d0][g1]

                            +PQ[c0]*delta[b0][g0]*delta[d0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_i * a_l * (
                            -PA_x*PB_0*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_i * (
                            +PA_x*PB_0*PQ[c0]*delta[d0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_0*QD_y*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_i * (
                            +QC_0*delta[d0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + S2 * inv_S1 * inv_S4 * a_i * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]
                        )

                    )

                    + F5_t[2] * (

                        S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -PA_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            +delta[b0][g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + inv_S4 * inv_S4 * a_i * a_l * (
                            +QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g1]*delta[c0][g0])

                            +delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g1]*delta[d0][g0]) + QD_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0])

                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PB_0*PQ[c0]*delta[d0][g1]

                            +PA_x*PB_0*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            -delta[d0][g1]*(PQ[c0] + QC_0)*(PA_x*PQ[b0] + PB_0*PQ[g0])

                            -PA_x*(PQ[c0]*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + PQ[d0]*(QC_0*delta[b0][g1] + QD_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0])) - PB_0*(PQ[c0]*(QD_0*delta[g0][g1] + QD_y*delta[d0][g0]) + PQ[d0]*(QC_0*delta[g0][g1] + QD_y*delta[c0][g0]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0]))

                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +QC_0*QD_y*(PQ[b0]*delta[d0][g0] + PQ[g0]*delta[b0][d0]) + QD_0*(PQ[b0]*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QD_y*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QD_y*delta[b0][c0]))

                            +PQ[b0]*PQ[g0]*(QC_0*delta[d0][g1] + QD_y*delta[c0][d0])

                            +delta[b0][g0]*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + (-1.0) * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*delta[c0][g0]*delta[d0][g1]

                            +delta[d0][g1]*(PQ[c0]*delta[b0][g0] + PQ[g0]*delta[b0][c0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PA_x*PB_0*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            -(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[c0]*delta[d0][g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[g0]*QC_0*QD_0*QD_y
                        )

                        + (-2.0) * S2 * S2 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[g0]*QC_0*delta[d0][g1]
                        )

                    )

                    + F5_t[3] * (

                        (-1.0) * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])

                            +PQ[c0]*(PA_x*(PQ[d0]*delta[b0][g1] + PQ[g1]*delta[b0][d0]) + PB_0*PQ[g1]*delta[d0][g0]) + PQ[d0]*PQ[g1]*(PA_x*delta[b0][c0] + PB_0*delta[c0][g0])

                            -PQ[c0]*PQ[d0]*PQ[g1]*delta[b0][g0]
                        )

                        + (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*(PQ[c0]*(QD_0*delta[g0][g1] + QD_y*delta[d0][g0]) + PQ[d0]*(QC_0*delta[g0][g1] + QD_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PQ[g0]*(PQ[c0]*(QD_0*delta[b0][g1] + QD_y*delta[b0][d0]) + PQ[d0]*(QC_0*delta[b0][g1] + QD_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]))

                            +PQ[b0]*PQ[g0]*(delta[c0][d0]*(PQ[g1] + QD_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PA_x*PB_0*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +(PA_x*PQ[b0] + PB_0*PQ[g0])*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PQ[b0]*PQ[g0]*(PQ[d0]*QC_0*QD_y + QD_0*(PQ[c0]*QD_y + PQ[g1]*QC_0))
                        )

                        + 2.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                            +PQ[b0]*PQ[c0]*PQ[g0]*delta[d0][g1]
                        )

                    )

                    + F5_t[4] * (

                        4.0 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            -PQ[c0]*PQ[d0]*PQ[g1]*(PA_x*PQ[b0] + PB_0*PQ[g0])
                        )

                        + 4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[g0]*(PQ[c0]*(PQ[d0]*QD_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                    )

                    + F5_t[5] * (

                        (-4.0) * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_l * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_il_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_il_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + atom_idx_l,
            hess_il_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + atom_idx_l * natoms + prim_cart_ao_to_atom_inds[i],
            hess_il_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeHessianSPPP_JK_0(double*         hess_xy,
                                double*         hess_yx,
                                const uint32_t  hess_cart_ind_0,
                                const uint32_t  hess_cart_ind_1,
                                const double    frac_exact_exchange,
                                const uint32_t* pair_inds_i_for_K_sp,
                                const uint32_t* pair_inds_k_for_K_sp,
                                const double*   D_ik_for_K_sp,
                                const uint32_t  pair_inds_count_for_K_sp,
                                const double*   s_prim_info,
                                const uint32_t* s_prim_aoinds,
                                const uint32_t  s_prim_count,
                                const double*   p_prim_info,
                                const uint32_t* p_prim_aoinds,
                                const uint32_t  p_prim_count,
                                const double    pp_max_D,
                                const double*   mat_D_full_AO,
                                const uint32_t  naos,
                                const double*   Q_K_sp,
                                const double*   Q_K_pp,
                                const uint32_t* D_inds_K_sp,
                                const uint32_t* D_inds_K_pp,
                                const uint32_t* pair_displs_K_sp,
                                const uint32_t* pair_displs_K_pp,
                                const uint32_t* pair_counts_K_sp,
                                const uint32_t* pair_counts_K_pp,
                                const double*   pair_data_K_sp,
                                const double*   pair_data_K_pp,
                                const uint32_t* prim_cart_ao_to_atom_inds,
                                const uint32_t  natoms,
                                const uint32_t* atom_ids_j,
                                const uint32_t* atom_ids_j_displs,
                                const uint32_t* atom_ids_j_counts,
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
    __shared__ uint32_t c0;
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ik = blockIdx.x;

    // we make sure that ik < pair_inds_count_for_K_sp when calling the kernel

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;

        i = pair_inds_i_for_K_sp[ik];
        k = pair_inds_k_for_K_sp[ik];

        count_i = pair_counts_K_sp[i];
        count_k = pair_counts_K_pp[k];

        displ_i = pair_displs_K_sp[i];
        displ_k = pair_displs_K_pp[k];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_k = p_prim_info[k / 3 + p_prim_count * 0];

        r_k[0] = p_prim_info[k / 3 + p_prim_count * 2];
        r_k[1] = p_prim_info[k / 3 + p_prim_count * 3];
        r_k[2] = p_prim_info[k / 3 + p_prim_count * 4];

        ik_factor_D = 2.0 * D_ik_for_K_sp[ik];

        c0 = k % 3;
    }

    __syncthreads();

    for (uint32_t atom_idx_j = 0; atom_idx_j < natoms; atom_idx_j++)
    {

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    __syncthreads();

    const auto atom_j_displ = atom_ids_j_displs[i * natoms + atom_idx_j];
    const auto atom_j_count = atom_ids_j_counts[i * natoms + atom_idx_j];

    for (uint32_t m = 0; m < (atom_j_count + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const auto j_raw = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PB_0, PB_x;
        uint32_t j_prim, j_cgto, b0;

        if (j_raw < atom_j_count)
        {
            // note non-coalesced read wrt j
            const auto j = atom_ids_j[displ_i + atom_j_displ + j_raw];

            Q_ij   = Q_K_sp[displ_i + j];

            j_prim = D_inds_K_sp[displ_i + j];

            j_cgto = p_prim_aoinds[(j_prim / 3) + p_prim_count * (j_prim % 3)];

            a_j = p_prim_info[j_prim / 3 + p_prim_count * 0];

            r_j[0] = p_prim_info[j_prim / 3 + p_prim_count * 2];
            r_j[1] = p_prim_info[j_prim / 3 + p_prim_count * 3];
            r_j[2] = p_prim_info[j_prim / 3 + p_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_sp[displ_i + j];

            PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);

            b0 = j_prim % 3;

            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

        }

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            // Q_kl == Q_K_pp[displ_k + l]
            if ((j_raw >= atom_j_count) || (l >= count_k) || (fabs(Q_ij * Q_K_pp[displ_k + l] * pp_max_D) <= eri_threshold))
            {
                break;
            }

            // const auto Q_kl = Q_K_pp[displ_k + l];

            const auto l_prim = D_inds_K_pp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_pp[displ_k + l];

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

            double d2 = 1.0;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, rho * d2 * r2_PQ, 5, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F5_t[1] *= d2;
                F5_t[2] *= d2 * d2;
                F5_t[3] *= d2 * d2 * d2;
                F5_t[4] *= d2 * d2 * d2 * d2;
                F5_t[5] *= d2 * d2 * d2 * d2 * d2;
            }

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

            const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

            // j-k Hessian

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    + F5_t[0] * (

                        inv_S1 * inv_S2 * a_j * a_k * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]

                            +delta[b0][g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + 2.0 * inv_S1 * a_j * a_k * (
                            +QC_0*QC_y*QD_0*delta[b0][g0]
                        )

                        + (-1.0) * inv_S1 * a_j * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S2 * a_j * a_k * (
                            +PB_0*PB_x*QC_0*delta[d0][g1]

                            +PB_0*PB_x*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + (-1.0) * inv_S2 * a_k * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]

                            +delta[b0][g0]*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + 4.0 * a_j * a_k * (
                            +PB_0*PB_x*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * a_j * (
                            +PB_0*PB_x*QD_0*delta[c0][g1]
                        )

                        + (-2.0) * a_k * (
                            +QC_0*QC_y*QD_0*delta[b0][g0]
                        )

                        + (
                            +QD_0*delta[b0][g0]*delta[c0][g1]
                        )

                    )

                    + F5_t[1] * (

                        (-1.0) * inv_S1 * inv_S4 * a_j * a_k * (
                            +QC_0*delta[b0][g0]*delta[d0][g1]

                            +delta[b0][g0]*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])
                        )

                        + inv_S2 * inv_S4 * a_j * a_k * (
                            +PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PB_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            -delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))
                        )

                        + (-2.0) * S1 * inv_S2 * inv_S4 * a_j * a_k * (
                            +PB_0*PB_x*delta[d0][g1]*(PQ[c0] + QC_0)

                            +PB_0*PB_x*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))
                        )

                        + S1 * inv_S2 * inv_S4 * a_k * (
                            +delta[b0][g0]*delta[c0][g1]*(PQ[d0] + QD_0)

                            +delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[d0][g1]*(PQ[c0] + QC_0))
                        )

                        + (-2.0) * S2 * inv_S1 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*QD_0*delta[b0][g0]
                        )

                        + S2 * inv_S1 * inv_S4 * a_j * (
                            +QD_0*delta[b0][g0]*delta[c0][g1]
                        )

                        + 2.0 * inv_S4 * a_j * a_k * (
                            +QC_0*delta[d0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])

                            +(PB_0*PQ[g0] + PB_x*PQ[b0])*(QC_y*delta[c0][d0] + QD_0*delta[c0][g1])

                            +QC_0*(PB_0*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PB_x*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1])) + QC_y*QD_0*(PB_0*delta[c0][g0] + PB_x*delta[b0][c0])

                            -delta[b0][g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + inv_S4 * a_j * (
                            -PB_0*delta[c0][g1]*delta[d0][g0]

                            -PB_x*delta[b0][d0]*delta[c0][g1]

                            +PQ[d0]*delta[b0][g0]*delta[c0][g1]
                        )

                        + 4.0 * S1 * inv_S4 * a_j * a_k * (
                            -PB_0*PB_x*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * inv_S4 * a_j * (
                            +PB_0*PB_x*PQ[d0]*delta[c0][g1]
                        )

                        + 2.0 * S1 * inv_S4 * a_k * (
                            +delta[b0][g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + (-1.0) * S1 * inv_S4 * (
                            +PQ[d0]*delta[b0][g0]*delta[c0][g1]
                        )

                        + 4.0 * S2 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*QD_0*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + (-2.0) * S2 * inv_S4 * a_j * (
                            +QD_0*delta[c0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                    )

                    + F5_t[2] * (

                        S1 * inv_S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PB_0*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            -PB_x*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            +delta[b0][g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + inv_S4 * inv_S4 * a_j * a_k * (
                            +QD_0*(delta[b0][c0]*delta[g0][g1] + delta[b0][g1]*delta[c0][g0])

                            +delta[b0][g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0) + delta[d0][g1]*(PQ[c0] + QC_0))

                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0])

                            +QC_0*(delta[b0][d0]*delta[g0][g1] + delta[b0][g1]*delta[d0][g0]) + QC_y*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0])
                        )

                        + 2.0 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_0*PB_x*PQ[c0]*delta[d0][g1]

                            +PB_0*PB_x*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + (-1.0) * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * (
                            +PQ[c0]*delta[b0][g0]*delta[d0][g1]

                            +delta[b0][g0]*(PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])
                        )

                        + 2.0 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            -delta[d0][g1]*(PQ[c0] + QC_0)*(PB_0*PQ[g0] + PB_x*PQ[b0])

                            -PB_0*(PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) - PB_x*(PQ[c0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[d0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]))

                            -(PB_0*PQ[g0] + PB_x*PQ[b0])*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +QC_0*QC_y*(PQ[b0]*delta[d0][g0] + PQ[g0]*delta[b0][d0]) + QD_0*(PQ[b0]*(PQ[g0]*delta[c0][g1] + QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]))

                            +delta[b0][g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))

                            +PQ[b0]*PQ[g0]*(QC_0*delta[d0][g1] + QC_y*delta[c0][d0])
                        )

                        + (-1.0) * S2 * inv_S4 * inv_S4 * a_j * (
                            +PQ[b0]*delta[c0][g1]*delta[d0][g0]

                            +delta[c0][g1]*(PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0])
                        )

                        + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_0*PB_x*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            -(PB_0*PQ[g0] + PB_x*PQ[b0])*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * (
                            +PQ[d0]*delta[c0][g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[g0]*QC_0*QC_y*QD_0
                        )

                        + (-2.0) * S2 * S2 * inv_S4 * inv_S4 * a_j * (
                            +PQ[b0]*PQ[g0]*QD_0*delta[c0][g1]
                        )

                    )

                    + F5_t[3] * (

                        (-1.0) * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*(delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0])

                            +PQ[c0]*(delta[b0][d0]*delta[g0][g1] + delta[b0][g0]*delta[d0][g1] + delta[b0][g1]*delta[d0][g0]) + PQ[d0]*(delta[b0][c0]*delta[g0][g1] + delta[b0][g0]*delta[c0][g1] + delta[b0][g1]*delta[c0][g0]) + PQ[g0]*(delta[b0][c0]*delta[d0][g1] + delta[b0][d0]*delta[c0][g1] + delta[b0][g1]*delta[c0][d0]) + PQ[g1]*(delta[b0][c0]*delta[d0][g0] + delta[b0][d0]*delta[c0][g0] + delta[b0][g0]*delta[c0][d0])
                        )

                        + 2.0 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PB_0*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +(PB_0*PQ[g0] + PB_x*PQ[b0])*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0])

                            +PB_x*PQ[c0]*(PQ[d0]*delta[b0][g1] + PQ[g1]*delta[b0][d0]) + PQ[g1]*(PB_0*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0]) + PB_x*PQ[d0]*delta[b0][c0])

                            -PQ[c0]*PQ[d0]*PQ[g1]*delta[b0][g0]
                        )

                        + (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*(PQ[c0]*(QC_y*delta[d0][g0] + QD_0*delta[g0][g1]) + PQ[d0]*(QC_0*delta[g0][g1] + QC_y*delta[c0][g0]) + PQ[g0]*delta[d0][g1]*(PQ[c0] + QC_0) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])) + PQ[g0]*(PQ[c0]*(QC_y*delta[b0][d0] + QD_0*delta[b0][g1]) + PQ[d0]*(QC_0*delta[b0][g1] + QC_y*delta[b0][c0]) + PQ[g1]*(QC_0*delta[b0][d0] + QD_0*delta[b0][c0]))

                            +PQ[b0]*PQ[g0]*(delta[c0][d0]*(PQ[g1] + QC_y) + delta[c0][g1]*(PQ[d0] + QD_0))

                            +delta[b0][g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PB_0*PB_x*PQ[c0]*PQ[d0]*PQ[g1]
                        )

                        + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +(PB_0*PQ[g0] + PB_x*PQ[b0])*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PQ[b0]*PQ[g0]*(PQ[g1]*QC_0*QD_0 + QC_y*(PQ[c0]*QD_0 + PQ[d0]*QC_0))
                        )

                        + 2.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                            +PQ[c0]*PQ[d0]*PQ[g1]*delta[b0][g0]
                        )

                        + 2.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * (
                            +PQ[b0]*PQ[d0]*PQ[g0]*delta[c0][g1]
                        )

                    )

                    + F5_t[4] * (

                        4.0 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            -PQ[c0]*PQ[d0]*PQ[g1]*(PB_0*PQ[g0] + PB_x*PQ[b0])
                        )

                        + 4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[g0]*(PQ[c0]*(PQ[d0]*QC_y + PQ[g1]*QD_0) + PQ[d0]*PQ[g1]*QC_0)
                        )

                        + 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*delta[g0][g1]

                            +PQ[g0]*(PQ[b0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[c0]*PQ[d0]*delta[b0][g1]) + PQ[g1]*(PQ[c0]*(PQ[b0]*delta[d0][g0] + PQ[d0]*delta[b0][g0] + PQ[g0]*delta[b0][d0]) + PQ[d0]*(PQ[b0]*delta[c0][g0] + PQ[g0]*delta[b0][c0]))
                        )

                    )

                    + F5_t[5] * (

                        (-4.0) * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_k * (
                            +PQ[b0]*PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                        )

                    )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
        }
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sp))
    {
        double hess_jk_xy = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                hess_jk_xy += ERIs[y][x];
            }
        }

        atomicAdd(
            hess_xy + atom_idx_j * natoms + prim_cart_ao_to_atom_inds[s_prim_count + k],
            hess_jk_xy * ik_factor_D * frac_exact_exchange);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[s_prim_count + k] * natoms + atom_idx_j,
            hess_jk_xy * ik_factor_D * frac_exact_exchange);
    }

    __syncthreads();

    }
}

}  // namespace gpu
