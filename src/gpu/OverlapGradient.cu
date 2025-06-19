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

#include "OverlapGradient.hpp"
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
computeOverlapGradientSS(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   ss_mat_W_local,
                         const uint32_t* first_inds_local,
                         const uint32_t* second_inds_local,
                         const uint32_t  ss_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds)
{
    // each thread computes a primitive S matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (ij < ss_prim_pair_count_local)
    {
        const auto i = first_inds_local[ij];
        const auto j = second_inds_local[ij];

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const double r_i[3] = {s_prim_info[i + s_prim_count * 2],
                               s_prim_info[i + s_prim_count * 3],
                               s_prim_info[i + s_prim_count * 4]};

        const auto a_j = s_prim_info[j + s_prim_count * 0];
        const auto c_j = s_prim_info[j + s_prim_count * 1];
        const double r_j[3] = {s_prim_info[j + s_prim_count * 2],
                               s_prim_info[j + s_prim_count * 3],
                               s_prim_info[j + s_prim_count * 4]};

        const auto r2_ij = (r_j[0] - r_i[0]) * (r_j[0] - r_i[0]) +
                           (r_j[1] - r_i[1]) * (r_j[1] - r_i[1]) +
                           (r_j[2] - r_i[2]) * (r_j[2] - r_i[2]);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto grad_i = S_ij_00 * (

            (

                2.0 * a_i * (
                    PA_x
                )

            )

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * ss_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[j], -grad_i * ss_mat_W_local[ij] * ij_factor);
    }
}

__global__ void
computeOverlapGradientSP(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   sp_mat_W_local,
                         const uint32_t* sp_first_inds_local,
                         const uint32_t* sp_second_inds_local,
                         const uint32_t  sp_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive S matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < sp_prim_pair_count_local)
    {
        const auto i = sp_first_inds_local[ij];
        const auto j = sp_second_inds_local[ij];

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const double r_i[3] = {s_prim_info[i + s_prim_count * 2],
                               s_prim_info[i + s_prim_count * 3],
                               s_prim_info[i + s_prim_count * 4]};

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const double r_j[3] = {p_prim_info[j / 3 + p_prim_count * 2],
                               p_prim_info[j / 3 + p_prim_count * 3],
                               p_prim_info[j / 3 + p_prim_count * 4]};

        const auto b0 = j % 3;

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = 2.0;

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        const auto grad_i = S_ij_00 * (

            (

                1.0 / (a_i + a_j) * a_i * (
                    delta[b0][grad_cart_ind]
                )

                + 2.0 * a_i * (
                    PA_x * PB_0
                )

            )

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * sp_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + j], -grad_i * sp_mat_W_local[ij] * ij_factor);
    }
}

__global__ void
computeOverlapGradientSD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   sd_mat_W_local,
                         const uint32_t* sd_first_inds_local,
                         const uint32_t* sd_second_inds_local,
                         const uint32_t  sd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  p_prim_count)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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
        const auto i = sd_first_inds_local[ij];
        const auto j = sd_second_inds_local[ij];

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const double r_i[3] = {s_prim_info[i + s_prim_count * 2],
                               s_prim_info[i + s_prim_count * 3],
                               s_prim_info[i + s_prim_count * 4]};

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const double r_j[3] = {d_prim_info[j / 6 + d_prim_count * 2],
                               d_prim_info[j / 6 + d_prim_count * 3],
                               d_prim_info[j / 6 + d_prim_count * 4]};

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = 2.0;

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const auto grad_i = S_ij_00 * (

            (

                1.0 / (a_i + a_j) * a_i * (
                    delta[b0][b1] * (PA_x)
                    + delta[b1][grad_cart_ind] * (PB_0)
                    + delta[b0][grad_cart_ind] * (PB_1)
                )

                + 2.0 * a_i * (
                    PA_x * PB_0 * PB_1
                )

            )

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * sd_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], -grad_i * sd_mat_W_local[ij] * ij_factor);
    }
}

__global__ void
computeOverlapGradientPP(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   pp_mat_W_local,
                         const uint32_t* pp_first_inds_local,
                         const uint32_t* pp_second_inds_local,
                         const uint32_t  pp_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive S matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < pp_prim_pair_count_local)
    {
        const auto i = pp_first_inds_local[ij];
        const auto j = pp_second_inds_local[ij];

        const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
        const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
        const double r_i[3] = {p_prim_info[i / 3 + p_prim_count * 2],
                               p_prim_info[i / 3 + p_prim_count * 3],
                               p_prim_info[i / 3 + p_prim_count * 4]};

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const double r_j[3] = {p_prim_info[j / 3 + p_prim_count * 2],
                               p_prim_info[j / 3 + p_prim_count * 3],
                               p_prim_info[j / 3 + p_prim_count * 4]};

        const auto a0 = i % 3;
        const auto b0 = j % 3;

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        const auto grad_i = S_ij_00 * (

            (

                (-1.0) * (
                    delta[a0][grad_cart_ind] * (PB_0)
                )

                + 1.0 / (a_i + a_j) * a_i * (
                    delta[b0][grad_cart_ind] * (PA_0)
                    + delta[a0][b0] * (PA_x)
                    + delta[a0][grad_cart_ind] * (PB_0)
                )

                + 2.0 * a_i * (
                    PA_0 * PA_x * PB_0
                )

            )

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + i], grad_i * pp_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + j], -grad_i * pp_mat_W_local[ij] * ij_factor);
    }
}

__global__ void
computeOverlapGradientPD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   pd_mat_W_local,
                         const uint32_t* pd_first_inds_local,
                         const uint32_t* pd_second_inds_local,
                         const uint32_t  pd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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
        const auto i = pd_first_inds_local[ij];
        const auto j = pd_second_inds_local[ij];

        const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
        const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
        const double r_i[3] = {p_prim_info[i / 3 + p_prim_count * 2],
                               p_prim_info[i / 3 + p_prim_count * 3],
                               p_prim_info[i / 3 + p_prim_count * 4]};

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const double r_j[3] = {d_prim_info[j / 6 + d_prim_count * 2],
                               d_prim_info[j / 6 + d_prim_count * 3],
                               d_prim_info[j / 6 + d_prim_count * 4]};

        const auto a0 = i % 3;

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = 2.0;

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const auto grad_i = S_ij_00 * (

            (

                (-0.5) / (a_i + a_j) * (
                    delta[a0][grad_cart_ind] * delta[b0][b1]
                )

                + (-1.0) * (
                    delta[a0][grad_cart_ind] * (PB_0 * PB_1)
                )

                + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                    (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                )

                + 1.0 / (a_i + a_j) * a_i * (
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

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + i], grad_i * pd_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], -grad_i * pd_mat_W_local[ij] * ij_factor);
    }
}

__global__ void
computeOverlapGradientDD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   dd_mat_W_local,
                         const uint32_t* dd_first_inds_local,
                         const uint32_t* dd_second_inds_local,
                         const uint32_t  dd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count,
                         const uint32_t  p_prim_count)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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
        const auto i = dd_first_inds_local[ij];
        const auto j = dd_second_inds_local[ij];

        const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
        const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
        const double r_i[3] = {d_prim_info[i / 6 + d_prim_count * 2],
                               d_prim_info[i / 6 + d_prim_count * 3],
                               d_prim_info[i / 6 + d_prim_count * 4]};

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const double r_j[3] = {d_prim_info[j / 6 + d_prim_count * 2],
                               d_prim_info[j / 6 + d_prim_count * 3],
                               d_prim_info[j / 6 + d_prim_count * 4]};

        const auto a0 = d_cart_inds[i % 6][0];
        const auto a1 = d_cart_inds[i % 6][1];

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const auto grad_i = S_ij_00 * (

            (

                (-0.5) / (a_i + a_j) * (
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

                + 1.0 / (a_i + a_j) * a_i * (
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

        );

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + i], grad_i * dd_mat_W_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], -grad_i * dd_mat_W_local[ij] * ij_factor);
    }
}

}  // namespace gpu
