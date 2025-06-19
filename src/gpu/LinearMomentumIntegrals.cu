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

#include "LinearMomentumIntegrals.hpp"
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
computeLinearMomentumSS(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const uint32_t* first_inds_local,
                        const uint32_t* second_inds_local,
                        const uint32_t  ss_prim_pair_count_local)
{
    // each thread computes a primitive matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (ij < ss_prim_pair_count_local)
    {
        const auto i = first_inds_local[ij];
        const auto j = second_inds_local[ij];

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = s_prim_info[j + s_prim_count * 0];
        const auto c_j = s_prim_info[j + s_prim_count * 1];
        const auto x_j = s_prim_info[j + s_prim_count * 2];
        const auto y_j = s_prim_info[j + s_prim_count * 3];
        const auto z_j = s_prim_info[j + s_prim_count * 4];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -2.0 * a_j * (PB_m)

                    );
        }
    }
}

__global__ void
computeLinearMomentumSP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const uint32_t* sp_first_inds_local,
                        const uint32_t* sp_second_inds_local,
                        const uint32_t  sp_prim_pair_count_local)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive matrix element

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
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const auto x_j = p_prim_info[j / 3 + p_prim_count * 2];
        const auto y_j = p_prim_info[j / 3 + p_prim_count * 3];
        const auto z_j = p_prim_info[j / 3 + p_prim_count * 4];

        const auto b0 = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -a_j / (a_i + a_j) * (delta[b0][m])

                    -2.0 * a_j * (PB_0 * PB_m)

                    + (delta[m][b0])

                    );
        }
    }
}

__global__ void
computeLinearMomentumSD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* sd_first_inds_local,
                        const uint32_t* sd_second_inds_local,
                        const uint32_t  sd_prim_pair_count_local)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive matrix element

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
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -a_j / (a_i + a_j) * (
                        delta[b1][m] * (PB_0)
                        + delta[b0][m] * (PB_1)
                        + delta[b0][b1] * (PB_m)
                    )

                    -2.0 * a_j * (
                        PB_0 * PB_1 * PB_m
                    )

                    + (
                        delta[m][b1] * (PB_0)
                        + delta[m][b0] * (PB_1)
                    )

                    );
        }
    }
}

__global__ void
computeLinearMomentumPP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const uint32_t* pp_first_inds_local,
                        const uint32_t* pp_second_inds_local,
                        const uint32_t  pp_prim_pair_count_local)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive matrix element

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
        const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
        const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
        const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const auto x_j = p_prim_info[j / 3 + p_prim_count * 2];
        const auto y_j = p_prim_info[j / 3 + p_prim_count * 3];
        const auto z_j = p_prim_info[j / 3 + p_prim_count * 4];

        const auto a0 = i % 3;
        const auto b0 = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -a_j / (a_i + a_j) * (
                        delta[b0][m] * (PA_0)
                        + delta[a0][m] * (PB_0)
                        + delta[a0][b0] * (PB_m)
                    )

                    -2.0 * a_j * (
                        PA_0 * PB_0 * PB_m
                    )

                    + (
                        delta[m][b0] * (PA_0)
                    )

                    );
        }
    }
}

__global__ void
computeLinearMomentumPD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* pd_first_inds_local,
                        const uint32_t* pd_second_inds_local,
                        const uint32_t  pd_prim_pair_count_local)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive matrix element

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
        const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
        const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
        const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];

        const auto a0 = i % 3;

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -0.5 * a_j / ( (a_i + a_j) * (a_i + a_j) ) * (
                        delta[a0][b0] * delta[b1][m]
                        + delta[a0][b1] * delta[b0][m]
                        + delta[a0][m] * delta[b0][b1]
                    )

                    -a_j / (a_i + a_j) * (
                        delta[b1][m] * (PA_0 * PB_0)
                        + delta[b0][m] * (PA_0 * PB_1)
                        + delta[b0][b1] * (PA_0 * PB_m)
                        + delta[a0][m] * (PB_0 * PB_1)
                        + delta[a0][b1] * (PB_0 * PB_m)
                        + delta[a0][b0] * (PB_1 * PB_m)
                    )

                    -2.0 * a_j * (
                        PA_0 * PB_0 * PB_1 * PB_m
                    )

                    + 0.5 / (a_i + a_j) * (
                        delta[a0][b0] * delta[m][b1]
                        + delta[a0][b1] * delta[m][b0]
                    )

                    + (
                        delta[m][b1] * (PA_0 * PB_0)
                        + delta[m][b0] * (PA_0 * PB_1)
                    )

                    );
        }
    }
}

__global__ void
computeLinearMomentumDD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* dd_first_inds_local,
                        const uint32_t* dd_second_inds_local,
                        const uint32_t  dd_prim_pair_count_local)
{
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive matrix element

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
        const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
        const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
        const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];

        const auto a0 = d_cart_inds[i % 6][0];
        const auto a1 = d_cart_inds[i % 6][1];

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto PB_m = (-a_i / (a_i + a_j)) * rij[m];

            mat_mu[m][ij] = S_ij_00 * (

                    -0.5 * a_j / ( (a_i + a_j) * (a_i + a_j) ) * (
                        delta[a1][b0] * delta[b1][m] * (PA_0)
                        + delta[a1][b1] * delta[b0][m] * (PA_0)
                        + delta[a1][m] * delta[b0][b1] * (PA_0)
                        + delta[a0][b0] * delta[b1][m] * (PA_1)
                        + delta[a0][b1] * delta[b0][m] * (PA_1)
                        + delta[a0][m] * delta[b0][b1] * (PA_1)
                        + delta[a0][a1] * delta[b1][m] * (PB_0)
                        + delta[a0][b1] * delta[a1][m] * (PB_0)
                        + delta[a0][m] * delta[a1][b1] * (PB_0)
                        + delta[a0][a1] * delta[b0][m] * (PB_1)
                        + delta[a0][b0] * delta[a1][m] * (PB_1)
                        + delta[a0][m] * delta[a1][b0] * (PB_1)
                        + delta[a0][a1] * delta[b0][b1] * (PB_m)
                        + delta[a0][b0] * delta[a1][b1] * (PB_m)
                        + delta[a0][b1] * delta[a1][b0] * (PB_m)
                    )

                    -a_j / (a_i + a_j) * (
                        delta[b1][m] * (PA_0 * PA_1 * PB_0)
                        + delta[b0][m] * (PA_0 * PA_1 * PB_1)
                        + delta[b0][b1] * (PA_0 * PA_1 * PB_m)
                        + delta[a1][m] * (PA_0 * PB_0 * PB_1)
                        + delta[a1][b1] * (PA_0 * PB_0 * PB_m)
                        + delta[a1][b0] * (PA_0 * PB_1 * PB_m)
                        + delta[a0][m] * (PA_1 * PB_0 * PB_1)
                        + delta[a0][b1] * (PA_1 * PB_0 * PB_m)
                        + delta[a0][b0] * (PA_1 * PB_1 * PB_m)
                        + delta[a0][a1] * (PB_0 * PB_1 * PB_m)
                    )

                    -2.0 * a_j * (
                        PA_0 * PA_1 * PB_0 * PB_1 * PB_m
                    )

                    + 0.5 / (a_i + a_j) * (
                        delta[a1][b0] * delta[m][b1] * (PA_0)
                        + delta[a1][b1] * delta[m][b0] * (PA_0)
                        + delta[a0][b0] * delta[m][b1] * (PA_1)
                        + delta[a0][b1] * delta[m][b0] * (PA_1)
                        + delta[a0][a1] * delta[m][b1] * (PB_0)
                        + delta[a0][a1] * delta[m][b0] * (PB_1)
                    )

                    + (
                        delta[m][b1] * (PA_0 * PA_1 * PB_0)
                        + delta[m][b0] * (PA_0 * PA_1 * PB_1)
                    )

                    );
        }
    }
}

}  // namespace gpu
