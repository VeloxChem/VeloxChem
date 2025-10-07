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
#include "EriExchangeGradientSSSP.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeGradientSSSP_I_0(double*         grad_x,
                                const int32_t  grad_cart_ind,
                                const double    frac_exact_exchange,
                                const int32_t* pair_inds_i_for_K_ss,
                                const int32_t* pair_inds_k_for_K_ss,
                                const double*   D_ik_for_K_ss,
                                const int32_t  pair_inds_count_for_K_ss,
                                const double*   s_prim_info,
                                const int32_t* s_prim_aoinds,
                                const int32_t  s_prim_count,
                                const double*   p_prim_info,
                                const int32_t* p_prim_aoinds,
                                const int32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const int32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_sp,
                                const int32_t* D_inds_K_ss,
                                const int32_t* D_inds_K_sp,
                                const int32_t* pair_displs_K_ss,
                                const int32_t* pair_displs_K_sp,
                                const int32_t* pair_counts_K_ss,
                                const int32_t* pair_counts_K_sp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_sp,
                                const int32_t* prim_cart_ao_to_atom_inds,
                                const double*   boys_func_table,
                                const double*   boys_func_ft,
                                const double    omega,
                                const double    eri_threshold)
{
    // each thread block scans over [i?|k?] and sum up to a primitive K matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM_Y_K][TILE_DIM_X_K + 1];
    __shared__ int32_t i, k, count_i, count_k, displ_i, displ_k;
    __shared__ double   a_i, r_i[3], a_k, r_k[3], ik_factor_D;
    __shared__ double   delta[3][3];

    const int32_t ik = blockIdx.x;

    double d2 = 1.0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        count_i = 0;
        count_k = 0;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        if (ik < pair_inds_count_for_K_ss)
        {
            i = rawValue(pair_inds_i_for_K_ss, ik);
            k = rawValue(pair_inds_k_for_K_ss, ik);

            count_i = rawValue(pair_counts_K_ss, i);
            count_k = rawValue(pair_counts_K_sp, k);

            displ_i = rawValue(pair_displs_K_ss, i);
            displ_k = rawValue(pair_displs_K_sp, k);

            a_i = rawValue(s_prim_info, i + s_prim_count * 0);

            r_i[0] = rawValue(s_prim_info, i + s_prim_count * 2);
            r_i[1] = rawValue(s_prim_info, i + s_prim_count * 3);
            r_i[2] = rawValue(s_prim_info, i + s_prim_count * 4);

            a_k = rawValue(s_prim_info, k + s_prim_count * 0);

            r_k[0] = rawValue(s_prim_info, k + s_prim_count * 2);
            r_k[1] = rawValue(s_prim_info, k + s_prim_count * 3);
            r_k[2] = rawValue(s_prim_info, k + s_prim_count * 4);

            ik_factor_D = (static_cast<double>(i != k) + 1.0) * rawValue(D_ik_for_K_ss, ik);


        }
    }

    __syncthreads();

    for (int32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const int32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        int32_t j_prim, j_cgto;

        if ((ik < pair_inds_count_for_K_ss) && (j < count_i))
        {
            Q_ij   = rawValue(Q_K_ss, displ_i + j);

            j_prim = rawValue(D_inds_K_ss, displ_i + j);

            j_cgto = rawValue(s_prim_aoinds, j_prim);

            a_j = rawValue(s_prim_info, j_prim + s_prim_count * 0);

            r_j[0] = rawValue(s_prim_info, j_prim + s_prim_count * 2);
            r_j[1] = rawValue(s_prim_info, j_prim + s_prim_count * 3);
            r_j[2] = rawValue(s_prim_info, j_prim + s_prim_count * 4);

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = rawValue(pair_data_K_ss, displ_i + j);

            PA_x = (a_j * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);



        }


        for (int32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const int32_t l = n * TILE_DIM_X_K + threadIdx.x;

            if ((ik >= pair_inds_count_for_K_ss) || (j >= count_i) || (l >= count_k) || (fabs(Q_ij * rawValue(Q_K_sp, displ_k + l) * sp_max_D) <= eri_threshold))
            {
                break;
            }

            const auto l_prim = rawValue(D_inds_K_sp, displ_k + l);

            const auto l_cgto = rawValue(p_prim_aoinds, (l_prim / 3) + p_prim_count * (l_prim % 3));

            const auto a_l = rawValue(p_prim_info, l_prim / 3 + p_prim_count * 0);

            const double r_l[3] = {rawValue(p_prim_info, l_prim / 3 + p_prim_count * 2),
                                   rawValue(p_prim_info, l_prim / 3 + p_prim_count * 3),
                                   rawValue(p_prim_info, l_prim / 3 + p_prim_count * 4)};

            const auto S_kl_00 = rawValue(pair_data_K_sp, displ_k + l);

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

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, rho * d2 * r2_PQ, 2, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F2_t[1] *= d2;
                F2_t[2] *= d2 * d2;
            }

            const double QC_x = (a_l * inv_S2) * (r_l[grad_cart_ind] - r_k[grad_cart_ind]);

            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);


            // mu grad

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                            F2_t[0] * 2.0 * a_i * (
                                + PA_x * QD_0
                            )

                            + F2_t[1] * (-2.0) * S1 * inv_S4 * a_i * (
                                PA_x * PQ[d0]
                            )

                            + F2_t[1] * 2.0 * S2 * inv_S4 * a_i * (
                                + PQ[grad_cart_ind] * QD_0
                            )

                            + F2_t[1] * inv_S4 * a_i * (
                                delta[d0][grad_cart_ind]
                            )

                            + F2_t[2] * (-2.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                                PQ[d0] * PQ[grad_cart_ind]
                            )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * rawValue(mat_D_full_AO, j_cgto * naos + l_cgto);

        }
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_ss))
    {
        double grad_i_x = 0.0;

        for (int32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (int32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                grad_i_x += ERIs[y][x];
            }
        }

        atomicAdd(grad_x + rawValue(prim_cart_ao_to_atom_inds, i), grad_i_x * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeGradientSSSP_K_0(double*         grad_x,
                                const int32_t  grad_cart_ind,
                                const double    frac_exact_exchange,
                                const int32_t* pair_inds_i_for_K_ss,
                                const int32_t* pair_inds_k_for_K_ss,
                                const double*   D_ik_for_K_ss,
                                const int32_t  pair_inds_count_for_K_ss,
                                const double*   s_prim_info,
                                const int32_t* s_prim_aoinds,
                                const int32_t  s_prim_count,
                                const double*   p_prim_info,
                                const int32_t* p_prim_aoinds,
                                const int32_t  p_prim_count,
                                const double    sp_max_D,
                                const double*   mat_D_full_AO,
                                const int32_t  naos,
                                const double*   Q_K_ss,
                                const double*   Q_K_sp,
                                const int32_t* D_inds_K_ss,
                                const int32_t* D_inds_K_sp,
                                const int32_t* pair_displs_K_ss,
                                const int32_t* pair_displs_K_sp,
                                const int32_t* pair_counts_K_ss,
                                const int32_t* pair_counts_K_sp,
                                const double*   pair_data_K_ss,
                                const double*   pair_data_K_sp,
                                const int32_t* prim_cart_ao_to_atom_inds,
                                const double*   boys_func_table,
                                const double*   boys_func_ft,
                                const double    omega,
                                const double    eri_threshold)
{
    // each thread block scans over [i?|k?] and sum up to a primitive K matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM_Y_K][TILE_DIM_X_K + 1];
    __shared__ int32_t i, k, count_i, count_k, displ_i, displ_k;
    __shared__ double   a_i, r_i[3], a_k, r_k[3], ik_factor_D;
    __shared__ double   delta[3][3];

    const int32_t ik = blockIdx.x;

    double d2 = 1.0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        count_i = 0;
        count_k = 0;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        if (ik < pair_inds_count_for_K_ss)
        {
            i = rawValue(pair_inds_i_for_K_ss, ik);
            k = rawValue(pair_inds_k_for_K_ss, ik);

            count_i = rawValue(pair_counts_K_ss, i);
            count_k = rawValue(pair_counts_K_sp, k);

            displ_i = rawValue(pair_displs_K_ss, i);
            displ_k = rawValue(pair_displs_K_sp, k);

            a_i = rawValue(s_prim_info, i + s_prim_count * 0);

            r_i[0] = rawValue(s_prim_info, i + s_prim_count * 2);
            r_i[1] = rawValue(s_prim_info, i + s_prim_count * 3);
            r_i[2] = rawValue(s_prim_info, i + s_prim_count * 4);

            a_k = rawValue(s_prim_info, k + s_prim_count * 0);

            r_k[0] = rawValue(s_prim_info, k + s_prim_count * 2);
            r_k[1] = rawValue(s_prim_info, k + s_prim_count * 3);
            r_k[2] = rawValue(s_prim_info, k + s_prim_count * 4);

            ik_factor_D = (static_cast<double>(i != k) + 1.0) * rawValue(D_ik_for_K_ss, ik);


        }
    }

    __syncthreads();

    for (int32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const int32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        int32_t j_prim, j_cgto;

        if ((ik < pair_inds_count_for_K_ss) && (j < count_i))
        {
            Q_ij   = rawValue(Q_K_ss, displ_i + j);

            j_prim = rawValue(D_inds_K_ss, displ_i + j);

            j_cgto = rawValue(s_prim_aoinds, j_prim);

            a_j = rawValue(s_prim_info, j_prim + s_prim_count * 0);

            r_j[0] = rawValue(s_prim_info, j_prim + s_prim_count * 2);
            r_j[1] = rawValue(s_prim_info, j_prim + s_prim_count * 3);
            r_j[2] = rawValue(s_prim_info, j_prim + s_prim_count * 4);

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = rawValue(pair_data_K_ss, displ_i + j);

            PA_x = (a_j * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);



        }


        for (int32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const int32_t l = n * TILE_DIM_X_K + threadIdx.x;

            if ((ik >= pair_inds_count_for_K_ss) || (j >= count_i) || (l >= count_k) || (fabs(Q_ij * rawValue(Q_K_sp, displ_k + l) * sp_max_D) <= eri_threshold))
            {
                break;
            }

            const auto l_prim = rawValue(D_inds_K_sp, displ_k + l);

            const auto l_cgto = rawValue(p_prim_aoinds, (l_prim / 3) + p_prim_count * (l_prim % 3));

            const auto a_l = rawValue(p_prim_info, l_prim / 3 + p_prim_count * 0);

            const double r_l[3] = {rawValue(p_prim_info, l_prim / 3 + p_prim_count * 2),
                                   rawValue(p_prim_info, l_prim / 3 + p_prim_count * 3),
                                   rawValue(p_prim_info, l_prim / 3 + p_prim_count * 4)};

            const auto S_kl_00 = rawValue(pair_data_K_sp, displ_k + l);

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

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, rho * d2 * r2_PQ, 2, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F2_t[1] *= d2;
                F2_t[2] *= d2 * d2;
            }

            const double QC_x = (a_l * inv_S2) * (r_l[grad_cart_ind] - r_k[grad_cart_ind]);

            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);


            // lambda grad

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                            F2_t[0] * 2.0 * a_k * (
                                + QC_x * QD_0
                            )

                            + F2_t[0] * inv_S2 * a_k * (
                                delta[d0][grad_cart_ind]
                            )

                            + F2_t[1] * (-1.0) * S1 * inv_S2 * inv_S4 * a_k * (
                                delta[d0][grad_cart_ind]
                            )

                            + F2_t[1] * (-2.0) * S1 * inv_S4 * a_k * (
                                PQ[d0] * QC_x
                                + PQ[grad_cart_ind] * QD_0
                            )

                            + F2_t[2] * 2.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                                PQ[d0] * PQ[grad_cart_ind]
                            )


                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * rawValue(mat_D_full_AO, l_cgto * naos + j_cgto);
        }
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_ss))
    {
        double grad_k_x = 0.0;

        for (int32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (int32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                grad_k_x += ERIs[y][x];
            }
        }

        atomicAdd(grad_x + rawValue(prim_cart_ao_to_atom_inds, k), grad_k_x * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}


}  // namespace gpu
