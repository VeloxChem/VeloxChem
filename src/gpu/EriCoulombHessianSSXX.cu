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
#include "EriCoulombHessianSSXX.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_II_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
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
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PA_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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

        // mu-mu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    (-2.0) * a_i * (
                        +delta[g0][g1]
                    )

                    + 2.0 * inv_S1 * a_i * a_i * (
                        +delta[g0][g1]
                    )

                    + 4.0 * a_i * a_i * (
                        +PA_x*PA_y
                    )

                )

                + F2_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                        +delta[g0][g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_i * (
                        +PA_x*PQ[g1]

                        +PA_y*PQ[g0]
                    )

                )

                + F2_t[2] * (

                    4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ii_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_JJ_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
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
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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

        // nu-nu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_j * a_j * (
                        +PB_x*PB_y
                    )

                    + (-2.0) * a_j * (
                        +delta[g0][g1]
                    )

                    + 2.0 * inv_S1 * a_j * a_j * (
                        +delta[g0][g1]
                    )

                )

                + F2_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_j * a_j * (
                        +delta[g0][g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_j * a_j * (
                        +PB_x*PQ[g1] + PB_y*PQ[g0]
                    )

                )

                + F2_t[2] * (

                    4.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_jj_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_jj_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[j], hess_jj_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_IJ_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* ss_first_inds,
                               const uint32_t* ss_second_inds,
                               const double*   ss_pair_data,
                               const uint32_t  ss_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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

        // mu-nu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_i * a_j * (
                        +PA_x*PB_y
                    )

                    + 2.0 * inv_S1 * a_i * a_j * (
                        +delta[g0][g1]
                    )

                )

                + F2_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                        +delta[g0][g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_j * (
                        +PA_x*PQ[g1] + PB_y*PQ[g0]
                    )

                )

                + F2_t[2] * (

                    4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ij_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[j],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[j] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_IK_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* ss_first_inds,
                               const uint32_t* ss_second_inds,
                               const double*   ss_pair_data,
                               const uint32_t  ss_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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


        const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // mu-lambda hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_i * a_k * (
                        +PA_x*QC_y
                    )

                )

                + F2_t[1] * (

                    (-4.0) * S1 * inv_S4 * a_i * a_k * (
                        +PA_x*PQ[g1]
                    )

                    + 2.0 * inv_S4 * a_i * a_k * (
                        +delta[g0][g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_k * (
                        +PQ[g0]*QC_y
                    )

                )

                + F2_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_k * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

        double hess_ik_xy = eri_ijkl * ss_mat_D[kl] * kl_factor;

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[k],
            hess_ik_xy * ij_factor_D * 2.0 * prefac_coulomb);

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ik_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ik_xy += ERIs[n][threadIdx.x];
        }

    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_JK_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* ss_first_inds,
                               const uint32_t* ss_second_inds,
                               const double*   ss_pair_data,
                               const uint32_t  ss_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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


        const auto QC_y = (a_l * inv_S2) * (r_l[g1] - r_k[g1]);

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // nu-lambda hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_j * a_k * (
                        +PB_x*QC_y
                    )

                )

                + F2_t[1] * (

                    4.0 * S1 * inv_S4 * a_j * a_k * (
                        -PB_x*PQ[g1]
                    )

                    + 2.0 * inv_S4 * a_j * a_k * (
                        +delta[g0][g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_j * a_k * (
                        +PQ[g0]*QC_y
                    )

                )

                + F2_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_j * a_k * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

        double hess_jk_xy = eri_ijkl * ss_mat_D[kl] * kl_factor;

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[j] * natoms + prim_cart_ao_to_atom_inds[k],
            hess_jk_xy * ij_factor_D * 2.0 * prefac_coulomb);

        /*
        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[k] * natoms + prim_cart_ao_to_atom_inds[j],
            hess_jk_xy * ij_factor_D * 2.0 * prefac_coulomb);
        */
    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_jk_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_jk_xy += ERIs[n][threadIdx.x];
        }

    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_IL_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* ss_first_inds,
                               const uint32_t* ss_second_inds,
                               const double*   ss_pair_data,
                               const uint32_t  ss_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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


        const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // mu-sigma hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_i * a_l * (
                        +PA_x*QD_y
                    )

                )

                + F2_t[1] * (

                    (-4.0) * S1 * inv_S4 * a_i * a_l * (
                        +PA_x*PQ[g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_l * (
                        +PQ[g0]*QD_y
                    )

                    + 2.0 * inv_S4 * a_i * a_l * (
                        +delta[g0][g1]
                    )

                )

                + F2_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_l * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

        double hess_il_xy = eri_ijkl * ss_mat_D[kl] * kl_factor;

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[l],
            hess_il_xy * ij_factor_D * 2.0 * prefac_coulomb);

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_il_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_il_xy += ERIs[n][threadIdx.x];
        }

    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSS_JL_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   ss_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   ss_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* ss_first_inds,
                               const uint32_t* ss_second_inds,
                               const double*   ss_pair_data,
                               const uint32_t  ss_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);



    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
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


        const auto QD_y = (-a_k * inv_S2) * (r_l[g1] - r_k[g1]);

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // nu-sigma hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F2_t[0] * (

                    4.0 * a_j * a_l * (
                        +PB_x*QD_y
                    )

                )

                + F2_t[1] * (

                    4.0 * S1 * inv_S4 * a_j * a_l * (
                        -PB_x*PQ[g1]
                    )

                    + 4.0 * S2 * inv_S4 * a_j * a_l * (
                        +PQ[g0]*QD_y
                    )

                    + 2.0 * inv_S4 * a_j * a_l * (
                        +delta[g0][g1]
                    )

                )

                + F2_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_j * a_l * (
                        +PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

        double hess_jl_xy = eri_ijkl * ss_mat_D[kl] * kl_factor;

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[j] * natoms + prim_cart_ao_to_atom_inds[l],
            hess_jl_xy * ij_factor_D * 2.0 * prefac_coulomb);

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_jl_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_jl_xy += ERIs[n][threadIdx.x];
        }

    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSP_II_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   sp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   sp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* sp_first_inds,
                               const uint32_t* sp_second_inds,
                               const double*   sp_pair_data,
                               const uint32_t  sp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PA_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= sp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * sp_mat_Q[kl] * sp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = sp_first_inds[kl];
        const auto l = sp_second_inds[kl];

        const auto a_k = s_prim_info[k + s_prim_count * 0];

        const double r_k[3] = {s_prim_info[k + s_prim_count * 2],
                               s_prim_info[k + s_prim_count * 3],
                               s_prim_info[k + s_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = sp_pair_data[kl];

        const auto d0 = l % 3;

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F3_t[4];

        gpu::computeBoysFunction(F3_t, S1 * S2 * inv_S4 * r2_PQ, 3, boys_func_table, boys_func_ft);

        const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

        double kl_factor = 2.0;

        // mu-mu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F3_t[0] * (

                    2.0 * inv_S1 * a_i * a_i * (
                        +QD_0*delta[g0][g1]
                    )

                    + 4.0 * a_i * a_i * (
                        +PA_x*PA_y*QD_0
                    )

                    + (-2.0) * a_i * (
                        +QD_0*delta[g0][g1]
                    )

                )

                + F3_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_i * (
                        +QD_0*delta[g0][g1]
                    )

                    + 2.0 * inv_S4 * a_i * a_i * (
                        -PQ[d0]*delta[g0][g1]

                        +PA_x*delta[d0][g1] + PA_y*delta[d0][g0]
                    )

                    + (-4.0) * S1 * inv_S4 * a_i * a_i * (
                        +PA_x*PA_y*PQ[d0]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_i * (
                        +QD_0*(PA_x*PQ[g1] + PA_y*PQ[g0])
                    )

                    + 2.0 * S1 * inv_S4 * a_i * (
                        +PQ[d0]*delta[g0][g1]
                    )

                )

                + F3_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                        +PA_x*PQ[d0]*PQ[g1]

                        +PA_y*PQ[d0]*PQ[g0]
                    )

                    + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                        +PQ[g0]*PQ[g1]*QD_0
                    )

                    + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_i * (
                        +PQ[d0]*delta[g0][g1]

                        +PQ[g0]*delta[d0][g1] + PQ[g1]*delta[d0][g0]
                    )

                )

                + F3_t[3] * (

                    (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_i * (
                        +PQ[d0]*PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * sp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ii_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSP_JJ_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   sp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   sp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* sp_first_inds,
                               const uint32_t* sp_second_inds,
                               const double*   sp_pair_data,
                               const uint32_t  sp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= sp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * sp_mat_Q[kl] * sp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = sp_first_inds[kl];
        const auto l = sp_second_inds[kl];

        const auto a_k = s_prim_info[k + s_prim_count * 0];

        const double r_k[3] = {s_prim_info[k + s_prim_count * 2],
                               s_prim_info[k + s_prim_count * 3],
                               s_prim_info[k + s_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = sp_pair_data[kl];

        const auto d0 = l % 3;

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F3_t[4];

        gpu::computeBoysFunction(F3_t, S1 * S2 * inv_S4 * r2_PQ, 3, boys_func_table, boys_func_ft);

        const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

        double kl_factor = 2.0;

        // nu-nu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F3_t[0] * (

                    2.0 * inv_S1 * a_j * a_j * (
                        +QD_0*delta[g0][g1]
                    )

                    + 4.0 * a_j * a_j * (
                        +PB_x*PB_y*QD_0
                    )

                    + (-2.0) * a_j * (
                        +QD_0*delta[g0][g1]
                    )

                )

                + F3_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_j * a_j * (
                        +QD_0*delta[g0][g1]
                    )

                    + 2.0 * inv_S4 * a_j * a_j * (
                        +PB_x*delta[d0][g1]

                        +PB_y*delta[d0][g0]

                        -PQ[d0]*delta[g0][g1]
                    )

                    + 4.0 * S1 * inv_S4 * a_j * a_j * (
                        -PB_x*PB_y*PQ[d0]
                    )

                    + 4.0 * S2 * inv_S4 * a_j * a_j * (
                        +QD_0*(PB_x*PQ[g1] + PB_y*PQ[g0])
                    )

                    + 2.0 * S1 * inv_S4 * a_j * (
                        +PQ[d0]*delta[g0][g1]
                    )

                )

                + F3_t[2] * (

                    4.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[d0]*(-PB_x*PQ[g1] - PB_y*PQ[g0])
                    )

                    + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[g0]*PQ[g1]*QD_0
                    )

                    + 2.0 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[d0]*delta[g0][g1]

                        +PQ[g0]*delta[d0][g1] + PQ[g1]*delta[d0][g0]
                    )

                )

                + F3_t[3] * (

                    (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[d0]*PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * sp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_jj_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_jj_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[j], hess_jj_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSSP_IJ_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   sp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   sp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* sp_first_inds,
                               const uint32_t* sp_second_inds,
                               const double*   sp_pair_data,
                               const uint32_t  sp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= sp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * sp_mat_Q[kl] * sp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = sp_first_inds[kl];
        const auto l = sp_second_inds[kl];

        const auto a_k = s_prim_info[k + s_prim_count * 0];

        const double r_k[3] = {s_prim_info[k + s_prim_count * 2],
                               s_prim_info[k + s_prim_count * 3],
                               s_prim_info[k + s_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = sp_pair_data[kl];

        const auto d0 = l % 3;

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F3_t[4];

        gpu::computeBoysFunction(F3_t, S1 * S2 * inv_S4 * r2_PQ, 3, boys_func_table, boys_func_ft);

        const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);

        double kl_factor = 2.0;

        // mu-nu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F3_t[0] * (

                    2.0 * inv_S1 * a_i * a_j * (
                        +QD_0*delta[g0][g1]
                    )

                    + 4.0 * a_i * a_j * (
                        +PA_x*PB_y*QD_0
                    )

                )

                + F3_t[1] * (

                    (-2.0) * S2 * inv_S1 * inv_S4 * a_i * a_j * (
                        +QD_0*delta[g0][g1]
                    )

                    + 2.0 * inv_S4 * a_i * a_j * (
                        +PB_y*delta[d0][g0]

                        -PQ[d0]*delta[g0][g1]

                        +PA_x*delta[d0][g1]
                    )

                    + 4.0 * S1 * inv_S4 * a_i * a_j * (
                        -PA_x*PB_y*PQ[d0]
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_j * (
                        +QD_0*(PA_x*PQ[g1] + PB_y*PQ[g0])
                    )

                )

                + F3_t[2] * (

                    (-4.0) * S1 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PA_x*PQ[d0]*PQ[g1]

                        +PB_y*PQ[d0]*PQ[g0]
                    )

                    + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PQ[g0]*PQ[g1]*QD_0
                    )

                    + 2.0 * S2 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PQ[d0]*delta[g0][g1]

                        +PQ[g0]*delta[d0][g1] + PQ[g1]*delta[d0][g0]
                    )

                )

                + F3_t[3] * (

                    (-4.0) * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PQ[d0]*PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * sp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ij_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[j],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[j] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSPP_II_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   pp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   pp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* pp_first_inds,
                               const uint32_t* pp_second_inds,
                               const double*   pp_pair_data,
                               const uint32_t  pp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PA_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PA_y = (a_j  * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= pp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * pp_mat_Q[kl] * pp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = pp_first_inds[kl];
        const auto l = pp_second_inds[kl];

        const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];

        const double r_k[3] = {p_prim_info[k / 3 + p_prim_count * 2],
                               p_prim_info[k / 3 + p_prim_count * 3],
                               p_prim_info[k / 3 + p_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = pp_pair_data[kl];

        const auto c0 = k % 3;
        const auto d0 = l % 3;

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

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // mu-mu hessian

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

                        +delta[g0][g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)

                        +delta[c0][d0]*(PA_x*PQ[g1] + PA_y*PQ[g0])
                    )

                    + 4.0 * S1 * inv_S4 * a_i * a_i * (
                        +PA_x*PA_y*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
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

                        -PA_y*PQ[c0]*delta[d0][g0] + PQ[d0]*(-PA_x*delta[c0][g1] - PA_y*delta[c0][g0])

                        +delta[c0][d0]*(-PA_x*PQ[g1] - PA_y*PQ[g0])

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
                        +PQ[g0]*PQ[g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
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

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * pp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ii_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ii_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[i], hess_ii_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSPP_JJ_0(double*         hess_xy,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   pp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   pp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* pp_first_inds,
                               const uint32_t* pp_second_inds,
                               const double*   pp_pair_data,
                               const uint32_t  pp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PB_x = (-a_i * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= pp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * pp_mat_Q[kl] * pp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = pp_first_inds[kl];
        const auto l = pp_second_inds[kl];

        const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];

        const double r_k[3] = {p_prim_info[k / 3 + p_prim_count * 2],
                               p_prim_info[k / 3 + p_prim_count * 3],
                               p_prim_info[k / 3 + p_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = pp_pair_data[kl];

        const auto c0 = k % 3;
        const auto d0 = l % 3;

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

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // nu-nu hessian

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                + F4_t[0] * (

                    2.0 * inv_S1 * a_j * a_j * (
                        +QC_0*QD_0*delta[g0][g1]
                    )

                    + 2.0 * inv_S2 * a_j * a_j * (
                        +PB_x*PB_y*delta[c0][d0]
                    )

                    + (-1.0) * inv_S2 * a_j * (
                        +delta[c0][d0]*delta[g0][g1]
                    )

                    + 4.0 * a_j * a_j * (
                        +PB_x*PB_y*QC_0*QD_0
                    )

                    + (-2.0) * a_j * (
                        +QC_0*QD_0*delta[g0][g1]
                    )

                    + inv_S1 * inv_S2 * a_j * a_j * (
                        +delta[c0][d0]*delta[g0][g1]
                    )

                )

                + F4_t[1] * (

                    (-1.0) * inv_S1 * inv_S4 * a_j * a_j * (
                        +delta[c0][d0]*delta[g0][g1]
                    )

                    + (-1.0) * inv_S2 * inv_S4 * a_j * a_j * (
                        +delta[c0][d0]*delta[g0][g1]
                    )

                    + (-2.0) * S1 * inv_S2 * inv_S4 * a_j * a_j * (
                        +PB_x*PB_y*delta[c0][d0]
                    )

                    + (-2.0) * S2 * inv_S1 * inv_S4 * a_j * a_j * (
                        +QC_0*QD_0*delta[g0][g1]
                    )

                    + 2.0 * inv_S4 * a_j * a_j * (
                        +delta[c0][d0]*(PB_x*PQ[g1] + PB_y*PQ[g0])

                        +PB_x*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PB_y*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                        +delta[g0][g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                    + 4.0 * S1 * inv_S4 * a_j * a_j * (
                        +PB_x*PB_y*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                    + 2.0 * S1 * inv_S4 * a_j * (
                        +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                    )

                    + 4.0 * S2 * inv_S4 * a_j * a_j * (
                        +QC_0*QD_0*(PB_x*PQ[g1] + PB_y*PQ[g0])
                    )

                    + S1 * inv_S2 * inv_S4 * a_j * (
                        +delta[c0][d0]*delta[g0][g1]
                    )

                )

                + F4_t[2] * (

                    2.0 * S1 * inv_S4 * inv_S4 * a_j * a_j * (
                        -PB_x*PQ[c0]*delta[d0][g1]

                        -PB_y*PQ[c0]*delta[d0][g0] + PQ[d0]*(-PB_x*delta[c0][g1] - PB_y*delta[c0][g0])

                        +delta[c0][d0]*(-PB_x*PQ[g1] - PB_y*PQ[g0])

                        +PQ[c0]*PQ[d0]*delta[g0][g1]
                    )

                    + 2.0 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +delta[g0][g1]*(PQ[c0]*QD_0 + PQ[d0]*QC_0)

                        +PQ[g0]*(QC_0*delta[d0][g1] + QD_0*delta[c0][g1]) + PQ[g1]*(QC_0*delta[d0][g0] + QD_0*delta[c0][g0])

                        +PQ[g0]*PQ[g1]*delta[c0][d0]
                    )

                    + 4.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PB_x*PB_y*PQ[c0]*PQ[d0]
                    )

                    + (-2.0) * S1 * S1 * inv_S4 * inv_S4 * a_j * (
                        +PQ[c0]*PQ[d0]*delta[g0][g1]
                    )

                    + 4.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        -(PB_x*PQ[g1] + PB_y*PQ[g0])*(PQ[c0]*QD_0 + PQ[d0]*QC_0)
                    )

                    + 4.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[g0]*PQ[g1]*QC_0*QD_0
                    )

                    + inv_S4 * inv_S4 * a_j * a_j * (
                        +delta[c0][d0]*delta[g0][g1] + delta[c0][g0]*delta[d0][g1] + delta[c0][g1]*delta[d0][g0]
                    )

                )

                + F4_t[3] * (

                    (-2.0) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[c0]*PQ[d0]*delta[g0][g1]

                        +PQ[g0]*(PQ[c0]*delta[d0][g1] + PQ[d0]*delta[c0][g1] + PQ[g1]*delta[c0][d0]) + PQ[g1]*(PQ[c0]*delta[d0][g0] + PQ[d0]*delta[c0][g0])
                    )

                    + 4.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[c0]*PQ[d0]*(PB_x*PQ[g1] + PB_y*PQ[g0])
                    )

                    + 4.0 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[g0]*PQ[g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                )

                + F4_t[4] * (

                    4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * a_j * (
                        +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * pp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_jj_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_jj_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(hess_xy + prim_cart_ao_to_atom_inds[j], hess_jj_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombHessianSSPP_IJ_0(double*         hess_xy,
                               double*         hess_yx,
                               const uint32_t  hess_cart_ind_0,
                               const uint32_t  hess_cart_ind_1,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   pp_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   pp_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* pp_first_inds,
                               const uint32_t* pp_second_inds,
                               const double*   pp_pair_data,
                               const uint32_t  pp_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const uint32_t  natoms,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];
    __shared__ uint32_t g0, g1;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PB_y;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        g0 = hess_cart_ind_0;
        g1 = hess_cart_ind_1;
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

        PA_x = (a_j  * inv_S1) * (r_j[g0] - r_i[g0]);
        PB_y = (-a_i * inv_S1) * (r_j[g1] - r_i[g1]);



    }

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= pp_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * pp_mat_Q[kl] * pp_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = pp_first_inds[kl];
        const auto l = pp_second_inds[kl];

        const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];

        const double r_k[3] = {p_prim_info[k / 3 + p_prim_count * 2],
                               p_prim_info[k / 3 + p_prim_count * 3],
                               p_prim_info[k / 3 + p_prim_count * 4]};

        const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];

        const double r_l[3] = {p_prim_info[l / 3 + p_prim_count * 2],
                               p_prim_info[l / 3 + p_prim_count * 3],
                               p_prim_info[l / 3 + p_prim_count * 4]};

        const auto S_kl_00 = pp_pair_data[kl];

        const auto c0 = k % 3;
        const auto d0 = l % 3;

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

        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // mu-nu hessian

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

                        +delta[g0][g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                    + 4.0 * S1 * inv_S4 * a_i * a_j * (
                        +PA_x*PB_y*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                    + 4.0 * S2 * inv_S4 * a_i * a_j * (
                        +QC_0*QD_0*(PA_x*PQ[g1] + PB_y*PQ[g0])
                    )

                )

                + F4_t[2] * (

                    2.0 * S1 * inv_S4 * inv_S4 * a_i * a_j * (
                        -PA_x*PQ[c0]*delta[d0][g1]

                        -PB_y*PQ[c0]*delta[d0][g0] + PQ[d0]*(-PA_x*delta[c0][g1] - PB_y*delta[c0][g0])

                        +delta[c0][d0]*(-PA_x*PQ[g1] - PB_y*PQ[g0])

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
                        +PQ[g0]*PQ[g1]*(-PQ[c0]*QD_0 - PQ[d0]*QC_0)
                    )

                )

                + F4_t[4] * (

                    4.0 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * a_j * (
                        +PQ[c0]*PQ[d0]*PQ[g0]*PQ[g1]
                    )

                )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * pp_mat_D[kl] * kl_factor;

    }

    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double hess_ij_xy = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            hess_ij_xy += ERIs[n][threadIdx.x];
        }

        atomicAdd(
            hess_xy + prim_cart_ao_to_atom_inds[i] * natoms + prim_cart_ao_to_atom_inds[j],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);

        atomicAdd(
            hess_yx + prim_cart_ao_to_atom_inds[j] * natoms + prim_cart_ao_to_atom_inds[i],
            hess_ij_xy * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

}  // namespace gpu
