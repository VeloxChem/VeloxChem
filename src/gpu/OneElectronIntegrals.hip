//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#include <hip/hip_runtime.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "BoysFuncTable.hpp"
#include "OneElectronIntegrals.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "StringFormat.hpp"

namespace gpu {  // gpu namespace

__device__ void
computeBoysFunction(double* values, const double fa, const uint32_t N, const double* bf_table, const double* ft)
{
    // Note: 847 = 121 * 7
    const double* bf_data = bf_table + N * 847;

    uint32_t pnt = (fa > 1.0e5) ? 1000000 : static_cast<uint32_t>(10.0 * fa + 0.5);

    if (pnt < 121)
    {
        const double w = fa - 0.1 * pnt;

        const double w2 = w * w;

        const double w4 = w2 * w2;

        values[N] = bf_data[pnt * 7 + 0] + bf_data[pnt * 7 + 1] * w + bf_data[pnt * 7 + 2] * w2 + bf_data[pnt * 7 + 3] * w2 * w

                    + bf_data[pnt * 7 + 4] * w4 + bf_data[pnt * 7 + 5] * w4 * w + bf_data[pnt * 7 + 6] * w4 * w2;

        const double f2a = fa + fa;

        const double fx = exp(-fa);

        for (uint32_t j = 0; j < N; j++)
        {
            values[N - j - 1] = ft[N - j - 1] * (f2a * values[N - j] + fx);
        }
    }
    else
    {
        const double fia = 1.0 / fa;

        double pf = 0.5 * fia;

        values[0] = MATH_CONST_HALF_SQRT_PI * sqrt(fia);

        if (pnt < 921)
        {
            const double fia2 = fia * fia;

            const double f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

            const double fx = exp(-fa);

            values[0] -= f * fx;

            const double rterm = pf * fx;

            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1] - rterm;

                pf += fia;
            }
        }
        else
        {
            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1];

                pf += fia;
            }
        }
    }
}

__global__ void
computeOverlapAndKineticEnergySS(double*         mat_S,
                                 double*         mat_T,
                                 const double*   s_prim_info,
                                 const uint32_t  s_prim_count,
                                 const uint32_t* first_inds_local,
                                 const uint32_t* second_inds_local,
                                 const uint32_t  ss_prim_pair_count_local)
{
    // each thread computes a primitive S/T matrix element

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00;

        mat_T[ij] = K_ij_00 * S_ij_00;
    }
}

__global__ void
computeOverlapAndKineticEnergySP(double*         mat_S,
                                 double*         mat_T,
                                 const double*   s_prim_info,
                                 const uint32_t  s_prim_count,
                                 const double*   p_prim_info,
                                 const uint32_t  p_prim_count,
                                 const uint32_t* sp_first_inds_local,
                                 const uint32_t* sp_second_inds_local,
                                 const uint32_t  sp_prim_pair_count_local)
{
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00 * PB_0;

        mat_T[ij] = S_ij_00 * (

                ksi * 2.0 * PB_0

                + K_ij_00 * PB_0

                );
    }
}

__global__ void
computeOverlapAndKineticEnergySD(double*         mat_S,
                                 double*         mat_T,
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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00 * (

                    0.5 / (a_i + a_j) * delta[b0][b1]

                    + PB_0 * PB_1

                );

        mat_T[ij] = S_ij_00 * (

                ksi * (

                    -1.0 / (a_i + a_j) * delta[b0][b1]

                    + 1.0 / (a_i + a_j) * delta[b0][b1]

                    + 4.0 * PB_0 * PB_1

                )

                + K_ij_00 * (

                    0.5 / (a_i + a_j) * delta[b0][b1]

                    + PB_0 * PB_1

                )

                );
    }
}

__global__ void
computeOverlapAndKineticEnergyPP(double*         mat_S,
                                 double*         mat_T,
                                 const double*   p_prim_info,
                                 const uint32_t  p_prim_count,
                                 const uint32_t* pp_first_inds_local,
                                 const uint32_t* pp_second_inds_local,
                                 const uint32_t  pp_prim_pair_count_local)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00 * (

                    0.5 / (a_i + a_j) * delta[b0][a0]

                    + PB_0 * PA_0

                );

        mat_T[ij] = S_ij_00 * (

                ksi * (

                    1.0 / (a_i + a_j) * delta[b0][a0]

                    + 4.0 * PB_0 * PA_0

                )

                + K_ij_00 * (

                    0.5 / (a_i + a_j) * delta[b0][a0]

                    + PB_0 * PA_0

                )

                );
    }
}

__global__ void
computeOverlapAndKineticEnergyPD(double*         mat_S,
                                 double*         mat_T,
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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00 * (

                    0.5 / (a_i + a_j) * (
                        delta[b1][a0] * (PB_0)
                        + delta[b0][a0] * (PB_1)
                        + delta[b0][b1] * (PA_0)
                    )

                    + (
                        PB_0 * PB_1 * PA_0
                    )

                );

        mat_T[ij] = S_ij_00 * (

                ksi * (

                    -1.0 / (a_i + a_j) * (
                        delta[b0][b1] * (PA_0)
                    )

                    + 2.0 / (a_i + a_j) * (
                        delta[b1][a0] * PB_0
                        + delta[b0][a0] * PB_1
                        + delta[b0][b1] * PA_0
                    )

                    + 6.0 * (
                        PB_0 * PB_1 * PA_0
                    )

                )

                + K_ij_00 * (

                    0.5 / (a_i + a_j) * (
                        delta[b1][a0] * (PB_0)
                        + delta[b0][a0] * (PB_1)
                        + delta[b0][b1] * (PA_0)
                    )

                    + (
                        PB_0 * PB_1 * PA_0
                    )

                )

                );
    }
}

__global__ void
computeOverlapAndKineticEnergyDD(double*         mat_S,
                                 double*         mat_T,
                                 const double*   d_prim_info,
                                 const uint32_t  d_prim_count,
                                 const uint32_t* dd_first_inds_local,
                                 const uint32_t* dd_second_inds_local,
                                 const uint32_t  dd_prim_pair_count_local)
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

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_S[ij] = S_ij_00 * (

                    0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                        delta[a0][a1] * delta[b0][b1]
                        + delta[b0][a0] * delta[b1][a1]
                        + delta[b0][a1] * delta[b1][a0]
                    )

                    + 0.5 / (a_i + a_j) * (
                        delta[a0][a1] * (PB_0 * PB_1)
                        + delta[b1][a1] * (PB_0 * PA_0)
                        + delta[b1][a0] * (PB_0 * PA_1)
                        + delta[b0][a1] * (PB_1 * PA_0)
                        + delta[b0][a0] * (PB_1 * PA_1)
                        + delta[b0][b1] * (PA_0 * PA_1)
                    )

                    + (
                        PB_0 * PB_1 * PA_0 * PA_1
                    )

                );

        mat_T[ij] = S_ij_00 * (

                ksi * (

                    1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (
                        delta[b0][a0] * delta[b1][a1]
                        + delta[b0][a1] * delta[b1][a0]
                    )

                    + 1.0 / (a_i + a_j) * (
                        delta[a0][a1] * (PB_0 * PB_1 * 2.0)
                        + delta[b1][a1] * (PB_0 * PA_0 * 3.0)
                        + delta[b1][a0] * (PB_0 * PA_1 * 3.0)
                        + delta[b0][a1] * (PB_1 * PA_0 * 3.0)
                        + delta[b0][a0] * (PB_1 * PA_1 * 3.0)
                        + delta[b0][b1] * (PA_0 * PA_1 * 2.0)
                    )

                    + 8.0 * (
                        PB_0 * PB_1 * PA_0 * PA_1
                    )

                )

                + K_ij_00 * (

                    0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                        delta[a0][a1] * delta[b0][b1]
                        + delta[b0][a0] * delta[b1][a1]
                        + delta[b0][a1] * delta[b1][a0]
                    )

                    + 0.5 / (a_i + a_j) * (
                        delta[a0][a1] * (PB_0 * PB_1)
                        + delta[b1][a1] * (PB_0 * PA_0)
                        + delta[b1][a0] * (PB_0 * PA_1)
                        + delta[b0][a1] * (PB_1 * PA_0)
                        + delta[b0][a0] * (PB_1 * PA_1)
                        + delta[b0][b1] * (PA_0 * PA_1)
                    )

                    + (
                        PB_0 * PB_1 * PA_0 * PA_1
                    )

                )

                );
    }
}

__global__ void
computeNuclearPotentialSS(double*         mat_V,
                          const double*   s_prim_info,
                          const uint32_t  s_prim_count,
                          const uint32_t* first_inds_local,
                          const uint32_t* second_inds_local,
                          const uint32_t  ss_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    // each thread computes a primitive S/T matrix element

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                  x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F0_t[1];

            gpu::computeBoysFunction(F0_t, (a_i + a_j) * r2_PC, 0, boys_func_table, boys_func_ft);

            V_ij += q_c * F0_t[0];
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeNuclearPotentialSP(double*         mat_V,
                          const double*   s_prim_info,
                          const uint32_t  s_prim_count,
                          const double*   p_prim_info,
                          const uint32_t  p_prim_count,
                          const uint32_t* sp_first_inds_local,
                          const uint32_t* sp_second_inds_local,
                          const uint32_t  sp_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                  y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F1_t[2];

            gpu::computeBoysFunction(F1_t, (a_i + a_j) * r2_PC, 1, boys_func_table, boys_func_ft);

            V_ij += q_c * (

                    F1_t[0] * PB_0

                    + F1_t[1] * PC[b0]

                    );
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeNuclearPotentialSD(double*         mat_V,
                          const double*   s_prim_info,
                          const uint32_t  s_prim_count,
                          const double*   d_prim_info,
                          const uint32_t  d_prim_count,
                          const uint32_t* sd_first_inds_local,
                          const uint32_t* sd_second_inds_local,
                          const uint32_t  sd_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                  y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * r2_PC, 2, boys_func_table, boys_func_ft);

            V_ij += q_c * (

                    F2_t[0] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][b1]
                        )

                        + (
                            PB_0 * PB_1
                        )

                    )

                    + F2_t[1] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b0][b1]
                        )

                        + (
                            PB_0 * PC[b1]
                            + PB_1 * PC[b0]
                        )

                    )

                    + F2_t[2] * (

                        (
                            PC[b0] * PC[b1]
                        )

                    )

                    );
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeNuclearPotentialPP(double*         mat_V,
                          const double*   p_prim_info,
                          const uint32_t  p_prim_count,
                          const uint32_t* pp_first_inds_local,
                          const uint32_t* pp_second_inds_local,
                          const uint32_t  pp_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                  y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * r2_PC, 2, boys_func_table, boys_func_ft);

            V_ij += q_c * (

                    F2_t[0] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][a0]
                        )

                        + (
                            PB_0 * PA_0
                        )

                    )

                    + F2_t[1] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b0][a0]
                        )

                        + (
                            PB_0 * PC[a0]
                            + PA_0 * PC[b0]
                        )

                    )

                    + F2_t[2] * (

                        (
                            PC[a0] * PC[b0]
                        )

                    )

                    );
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeNuclearPotentialPD(double*         mat_V,
                          const double*   p_prim_info,
                          const uint32_t  p_prim_count,
                          const double*   d_prim_info,
                          const uint32_t  d_prim_count,
                          const uint32_t* pd_first_inds_local,
                          const uint32_t* pd_second_inds_local,
                          const uint32_t  pd_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
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

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                  y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, (a_i + a_j) * r2_PC, 3, boys_func_table, boys_func_ft);

            V_ij += q_c * (

                    F3_t[0] * (

                        0.5 / (a_i + a_j) * (
                            delta[b1][a0] * (PB_0)
                            + delta[b0][a0] * (PB_1)
                            + delta[b0][b1] * (PA_0)
                        )

                        + (
                            PB_0 * PB_1 * PA_0
                        )

                    )

                    + F3_t[1] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b1][a0] * (PB_0)
                            + delta[b0][a0] * (PB_1)
                            + delta[b0][b1] * (PA_0)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PC[a0])
                            + delta[b1][a0] * (PC[b0])
                            + delta[b0][a0] * (PC[b1])
                        )

                        + (
                            PB_0 * PB_1 * PC[a0]
                            + PB_0 * PA_0 * PC[b1]
                            + PB_1 * PA_0 * PC[b0]
                        )

                    )

                    + F3_t[2] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PC[a0])
                            + delta[b1][a0] * (PC[b0])
                            + delta[b0][a0] * (PC[b1])
                        )

                        + (
                            PB_0 * PC[a0] * PC[b1]
                            + PB_1 * PC[a0] * PC[b0]
                            + PA_0 * PC[b0] * PC[b1]
                        )

                    )

                    + F3_t[3] * (

                        (
                            PC[a0] * PC[b0] * PC[b1]
                        )

                    )

                    );
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeNuclearPotentialDD(double*         mat_V,
                          const double*   d_prim_info,
                          const uint32_t  d_prim_count,
                          const uint32_t* dd_first_inds_local,
                          const uint32_t* dd_second_inds_local,
                          const uint32_t  dd_prim_pair_count_local,
                          const double*   points_info,
                          const uint32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
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

        // TODO: improve memory access pattern

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

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        double V_ij = 0.0;

        // TODO: use 2D block to scan the points

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                  y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                  z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F4_t[5];

            gpu::computeBoysFunction(F4_t, (a_i + a_j) * r2_PC, 4, boys_func_table, boys_func_ft);

            V_ij += q_c * (

                    F4_t[0] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[b0][b1]
                            + delta[b0][a0] * delta[b1][a1]
                            + delta[b0][a1] * delta[b1][a0]
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[a0][a1] * (PB_0 * PB_1)
                            + delta[b1][a1] * (PB_0 * PA_0)
                            + delta[b1][a0] * (PB_0 * PA_1)
                            + delta[b0][a1] * (PB_1 * PA_0)
                            + delta[b0][a0] * (PB_1 * PA_1)
                            + delta[b0][b1] * (PA_0 * PA_1)
                        )

                        + (
                            PB_0 * PB_1 * PA_0 * PA_1
                        )

                    )

                    + F4_t[1] * (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[b0][b1] * 2.0
                            + delta[b0][a0] * delta[b1][a1] * 2.0
                            + delta[b0][a1] * delta[b1][a0] * 2.0
                        )

                        -0.5 / (a_i + a_j) * (
                            delta[a0][a1] * (PB_0 * PB_1)
                            + delta[b1][a1] * (PB_0 * PA_0)
                            + delta[b1][a0] * (PB_0 * PA_1)
                            + delta[b0][a1] * (PB_1 * PA_0)
                            + delta[b0][a0] * (PB_1 * PA_1)
                            + delta[b0][b1] * (PA_0 * PA_1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][a1] * (PB_0 * PC[a0] + PA_0 * PC[b0])
                            + delta[b1][a0] * (PB_0 * PC[a1] + PA_1 * PC[b0])
                            + delta[a0][a1] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                            + delta[b0][a1] * (PB_1 * PC[a0] + PA_0 * PC[b1])
                            + delta[b0][a0] * (PB_1 * PC[a1] + PA_1 * PC[b1])
                            + delta[b0][b1] * (PA_0 * PC[a1] + PA_1 * PC[a0])
                        )

                        + (
                            PB_0 * PB_1 * PA_0 * PC[a1]
                            + PB_0 * PB_1 * PA_1 * PC[a0]
                            + PB_0 * PA_0 * PA_1 * PC[b1]
                            + PB_1 * PA_0 * PA_1 * PC[b0]
                        )

                    )

                    + F4_t[2] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b1][a1] * (PB_0 * PC[a0] + PA_0 * PC[b0])
                            + delta[b1][a0] * (PB_0 * PC[a1] + PA_1 * PC[b0])
                            + delta[a0][a1] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                            + delta[b0][a1] * (PB_1 * PC[a0] + PA_0 * PC[b1])
                            + delta[b0][a0] * (PB_1 * PC[a1] + PA_1 * PC[b1])
                            + delta[b0][b1] * (PA_0 * PC[a1] + PA_1 * PC[a0])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[b0][b1]
                            + delta[b0][a0] * delta[b1][a1]
                            + delta[b0][a1] * delta[b1][a0]
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PC[a0] * PC[a1])
                            + delta[b1][a1] * (PC[a0] * PC[b0])
                            + delta[b0][a1] * (PC[a0] * PC[b1])
                            + delta[b1][a0] * (PC[a1] * PC[b0])
                            + delta[b0][a0] * (PC[a1] * PC[b1])
                            + delta[a0][a1] * (PC[b0] * PC[b1])
                        )

                        + (
                            PB_0 * PB_1 * PC[a0] * PC[a1]
                            + PB_0 * PA_0 * PC[a1] * PC[b1]
                            + PB_0 * PA_1 * PC[a0] * PC[b1]
                            + PB_1 * PA_0 * PC[a1] * PC[b0]
                            + PB_1 * PA_1 * PC[a0] * PC[b0]
                            + PA_0 * PA_1 * PC[b0] * PC[b1]
                        )

                    )

                    + F4_t[3] * (

                        -0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PC[a0] * PC[a1])
                            + delta[b1][a1] * (PC[a0] * PC[b0])
                            + delta[b0][a1] * (PC[a0] * PC[b1])
                            + delta[b1][a0] * (PC[a1] * PC[b0])
                            + delta[b0][a0] * (PC[a1] * PC[b1])
                            + delta[a0][a1] * (PC[b0] * PC[b1])
                        )

                        + (
                            PB_0 * PC[a0] * PC[a1] * PC[b1]
                            + PB_1 * PC[a0] * PC[a1] * PC[b0]
                            + PA_0 * PC[a1] * PC[b0] * PC[b1]
                            + PA_1 * PC[a0] * PC[b0] * PC[b1]
                        )

                    )

                    + F4_t[4] * (

                        (
                            PC[a0] * PC[a1] * PC[b0] * PC[b1]
                        )

                    )

                    );
        }

        mat_V[ij] = theta * S_ij_00 * V_ij;
    }
}

__global__ void
computeQMatrixSS(double*         mat_Q,
                 const double*   s_prim_info,
                 const uint32_t  s_prim_count,
                 const uint32_t* ss_first_inds_local,
                 const uint32_t* ss_second_inds_local,
                 const uint32_t  ss_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (ij < ss_prim_pair_count_local)
    {
        const auto i = ss_first_inds_local[ij];
        const auto j = ss_second_inds_local[ij];

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F0_t[1];

        gpu::computeBoysFunction(F0_t, 0.0, 0, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * F0_t[0];

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

__global__ void
computeQMatrixSP(double*         mat_Q,
                 const double*   s_prim_info,
                 const uint32_t  s_prim_count,
                 const double*   p_prim_info,
                 const uint32_t  p_prim_count,
                 const uint32_t* sp_first_inds_local,
                 const uint32_t* sp_second_inds_local,
                 const uint32_t  sp_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F1_t[2];

        gpu::computeBoysFunction(F1_t, 0.0, 1, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto PB_0 = (-a_i / S1) * rij[b0];

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * (

                    F1_t[0] * PB_0 * PB_0

                    + F1_t[1] * 0.25 / S1

                );

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

__global__ void
computeQMatrixSD(double*         mat_Q,
                 const double*   s_prim_info,
                 const uint32_t  s_prim_count,
                 const double*   d_prim_info,
                 const uint32_t  d_prim_count,
                 const uint32_t* sd_first_inds_local,
                 const uint32_t* sd_second_inds_local,
                 const uint32_t  sd_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F2_t[3];

        gpu::computeBoysFunction(F2_t, 0.0, 2, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto PB_0 = (-a_i / S1) * rij[b0];
        const auto PB_1 = (-a_i / S1) * rij[b1];

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * (

                    F2_t[0] * (

                        0.25 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1]
                        )

                        + 0.5 / S1 * (
                            delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                        )

                        + (
                            PB_0 * PB_0 * PB_1 * PB_1
                        )

                    )

                    + F2_t[1] * (

                        0.125 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1] * (1.0 * (-2.0))
                        )

                        + 0.25 / S1 * (
                            + 1.0 * (PB_0 * PB_0)
                            + 1.0 * (PB_1 * PB_1)
                        )

                    )

                    + F2_t[2] * (

                        0.0625 / ( S1 * S1 ) * (
                            1.0
                            + delta[b0][b1] * delta[b0][b1] * 2.0
                        )

                    )

                );

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

__global__ void
computeQMatrixPP(double*         mat_Q,
                 const double*   p_prim_info,
                 const uint32_t  p_prim_count,
                 const uint32_t* pp_first_inds_local,
                 const uint32_t* pp_second_inds_local,
                 const uint32_t  pp_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    __shared__ double   delta[3][3];

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F2_t[3];

        gpu::computeBoysFunction(F2_t, 0.0, 2, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto PA_0 = (a_j / S1) * rij[a0];
        const auto PB_0 = (-a_i / S1) * rij[b0];

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * (

                    F2_t[0] * (

                        0.5 / S1 * (
                            delta[a0][b0] * ( + PB_0 * PA_0 * 2.0)
                        )

                        + (

                            + PB_0 * PB_0 * PA_0 * PA_0
                        )

                        + 0.25 / ( S1 * S1 ) * (
                            delta[a0][b0] * delta[a0][b0]
                        )

                    )

                    + F2_t[1] * (

                        0.125 / ( S1 * S1 ) * (
                            delta[a0][b0] * delta[a0][b0] * (1.0 * (-2.0))
                        )

                        + 0.25 / S1 * (
                            1.0 * ( + PA_0 * PA_0)
                            + 1.0 * (PB_0 * PB_0)
                        )

                    )

                    + F2_t[2] * (

                        0.0625 / ( S1 * S1 ) * (
                            1.0
                            + delta[a0][b0] * delta[a0][b0] * 2.0
                        )

                    )

                );

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

__global__ void
computeQMatrixPD(double*         mat_Q,
                 const double*   p_prim_info,
                 const uint32_t  p_prim_count,
                 const double*   d_prim_info,
                 const uint32_t  d_prim_count,
                 const uint32_t* pd_first_inds_local,
                 const uint32_t* pd_second_inds_local,
                 const uint32_t  pd_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F3_t[4];

        gpu::computeBoysFunction(F3_t, 0.0, 3, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto PA_0 = (a_j / S1) * rij[a0];
        const auto PB_0 = (-a_i / S1) * rij[b0];
        const auto PB_1 = (-a_i / S1) * rij[b1];

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * (

                    F3_t[0] * (

                        0.25 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0)
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                            + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                        )

                        + 0.5 / S1 * (
                            delta[b0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                            + delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                            + delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                        )

                        + (

                            + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0
                        )

                    )

                    + F3_t[1] * (

                        0.125 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * (-2.0))
                            + delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][b0] * ( + PB_0 * PA_0 * 2.0)
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * (-2.0))
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                            + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-1.0))
                            + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                            + 1.0 * delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                            + delta[a0][b1] * delta[a0][b0] * (PB_0 * PB_1)
                        )

                        + 0.25 / S1 * (
                            1.0 * ( + PB_0 * PB_0 * PA_0 * PA_0)
                            + 1.0 * ( + PB_1 * PB_1 * PA_0 * PA_0)
                            + 1.0 * (PB_0 * PB_0 * PB_1 * PB_1)
                        )

                        + 0.0625 / ( S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1]
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 4.0
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1]
                            + delta[a0][b0] * delta[a0][b0]
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                            + 1.0 * delta[a0][b1] * delta[a0][b1]
                        )

                    )

                    + F3_t[2] * (

                        0.03125 / ( S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (1.0 * (-8.0))
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[a0][b0] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * (1.0 * (-2.0))
                        )

                        + 0.0625 / ( S1 * S1 ) * (
                            1.0 * ( + PA_0 * PA_0)
                            + delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * (-1.0))
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0)
                            + delta[a0][b1] * ( + PB_1 * PA_0)
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * (-1.0))
                            + 1.0 * (PB_0 * PB_0)
                            + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * 2.0)
                            + 1.0 * (PB_1 * PB_1)
                            + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * 2.0)
                        )

                    )

                    + F3_t[3] * (

                        0.015625 / ( S1 * S1 * S1 ) * (
                            1.0
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * 2.0
                            + delta[a0][b0] * delta[a0][b0] * 2.0
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 5.0
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                            + delta[a0][b1] * delta[a0][b1]
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                            + 1.0 * delta[a0][b1] * delta[a0][b1]
                        )

                    )

                );

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

__global__ void
computeQMatrixDD(double*         mat_Q,
                 const double*   d_prim_info,
                 const uint32_t  d_prim_count,
                 const uint32_t* dd_first_inds_local,
                 const uint32_t* dd_second_inds_local,
                 const uint32_t  dd_prim_pair_count_local,
                 const double*   boys_func_table,
                 const double*   boys_func_ft)
{
    // each thread computes a primitive Q matrix element

    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

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

        const auto S1 = (a_i + a_j);

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto Lambda = sqrt(2.0 * S1 / MATH_CONST_PI);

        double F4_t[5];

        gpu::computeBoysFunction(F4_t, 0.0, 4, boys_func_table, boys_func_ft);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / S1, 1.5) * exp(-a_i * a_j / S1 * r2_ij);

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto PA_0 = (a_j / S1) * rij[a0];
        const auto PA_1 = (a_j / S1) * rij[a1];
        const auto PB_0 = (-a_i / S1) * rij[b0];
        const auto PB_1 = (-a_i / S1) * rij[b1];

        const double eri_ijij = Lambda * S_ij_00 * S_ij_00 * (

                    F4_t[0] * (

                        0.125 / ( S1 * S1 * S1 ) * (
                            delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * 2.0)
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                        )

                        + 0.25 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * PA_1 * PA_1)
                            + delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PB_0 * PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b1] * ( + PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PB_1 * PA_0 * PA_0)
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PB_0 * PA_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[b0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                            + delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                            + delta[a1][b0] * delta[a0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                            + delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a1][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                            + delta[a0][a1] * delta[a1][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PB_0 * PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b1] * ( + PB_0 * PB_1 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PB_1 * PA_1 * PA_1)
                            + delta[a0][a1] * delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1)
                        )

                        + 0.5 / S1 * (
                            delta[b0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a1][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][a1] * ( + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                        )

                        + (

                            + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1
                        )

                        + 0.0625 / ( S1 * S1 * S1 * S1 ) * (
                            delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1]
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 2.0
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 2.0
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1]
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 2.0
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1]
                        )

                    )

                    + F4_t[1] * (

                        0.03125 / ( S1 * S1 * S1 * S1 ) * (
                            delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * (1.0 * (-4.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * (1.0 * (-8.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * (1.0 * (-8.0))
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (1.0 * (-4.0))
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * (1.0 * (-8.0))
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * (1.0 * (-4.0))
                        )

                        + 0.0625 / ( S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0)
                            + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_0 * 4.0)
                            + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * ( + PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b0] * ( + PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * ( + PA_0 * PA_0)
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * ( + PA_0 * PA_0)
                            + delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_1 * (-4.0))
                            + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_1 * (-2.0))
                            + delta[a0][b1] * delta[a1][b0] * delta[b0][b1] * ( + PA_0 * PA_1)
                            + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_0 * PA_1 * (-2.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_0 * PA_1)
                            + delta[a0][b0] * delta[a1][b0] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * ( + PA_0 * PA_1)
                            + delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * ( + PA_0 * PA_1)
                            + 1.0 * delta[a0][b1] * delta[a1][b1] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * (-2.0))
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * (-4.0))
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b1] * ( + PB_0 * PA_0)
                            + delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][a1] * delta[a1][b0] * ( + PB_0 * PA_0 * 2.0)
                            + 1.0 * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * (-2.0))
                            + delta[a1][b1] * delta[a1][b0] * delta[a0][b1] * ( + PB_0 * PA_0)
                            + delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * ( + PB_1 * PA_0 * (-2.0))
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PA_0)
                            + delta[a0][a1] * delta[a1][b1] * ( + PB_1 * PA_0 * 2.0)
                            + 1.0 * delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a0][b0] * delta[a1][b1] * ( + PB_1 * PA_0)
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PA_0)
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_0 * (-5.0))
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_1 * PA_1)
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_1 * PA_1 * 4.0)
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b0] * ( + PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * ( + PA_1 * PA_1)
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * ( + PA_1 * PA_1)
                            + 1.0 * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_1 * (-2.0))
                            + delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PA_1 * (-2.0))
                            + delta[a0][b1] * delta[a0][a1] * delta[b0][b1] * ( + PB_0 * PA_1)
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PA_1)
                            + delta[a0][a1] * delta[a0][b0] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_1)
                            + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PA_1 * (-4.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b1] * ( + PB_0 * PA_1)
                            + 1.0 * delta[a1][b0] * delta[b0][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * ( + PB_1 * PA_1)
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_1 * (-2.0))
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PA_1)
                            + delta[a0][a1] * delta[a0][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * ( + PB_1 * PA_1)
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * ( + PB_1 * PA_1 * (-5.0))
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * (-4.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-2.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-1.0))
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0)
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * 4.0)
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0)
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0)
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0)
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                            + 1.0 * delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * 2.0)
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_1)
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * (PB_0 * PB_1)
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1)
                            + 1.0 * delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                            + 1.0 * delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1)
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * 4.0)
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1)
                            + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1)
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1)
                            + 1.0 * delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                        )

                        + 0.125 / ( S1 * S1 ) * (
                            delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * PA_1 * PA_1 * (-2.0))
                            + delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * ( + PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * PA_0 * PA_1 * (-2.0))
                            + 1.0 * delta[a1][b1] * ( + PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PB_0 * PA_0 * PA_0 * (-2.0))
                            + 1.0 * delta[b0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a1][b1] * ( + PB_0 * PB_1 * PA_0 * PA_0 * (-1.0))
                            + delta[a1][b1] * delta[a1][b0] * ( + PB_0 * PB_1 * PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PB_1 * PA_0 * PA_0 * (-2.0))
                            + delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * ( + PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * PA_1 * PA_1 * (-2.0))
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b1] * delta[a1][b1] * ( + PB_0 * PB_0 * PA_0 * PA_1 * (-1.0))
                            + delta[a0][a1] * ( + PB_0 * PB_0 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b1] * delta[a0][b1] * ( + PB_0 * PB_0 * PA_0 * PA_1)
                            + delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1 * (-2.0))
                            + delta[a0][b1] * delta[a1][b0] * ( + PB_0 * PB_1 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1 * (-2.0))
                            + delta[a1][b1] * delta[a0][b0] * ( + PB_0 * PB_1 * PA_0 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PB_1 * PA_0 * PA_1 * (-1.0))
                            + delta[a0][a1] * ( + PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b0] * ( + PB_1 * PB_1 * PA_0 * PA_1)
                            + 1.0 * delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                            + 1.0 * delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PB_0 * PA_1 * PA_1 * (-2.0))
                            + 1.0 * delta[b0][b1] * ( + PB_0 * PB_1 * PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b1] * ( + PB_0 * PB_1 * PA_1 * PA_1 * (-1.0))
                            + delta[a0][b1] * delta[a0][b0] * ( + PB_0 * PB_1 * PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PB_1 * PA_1 * PA_1 * (-2.0))
                            + 1.0 * delta[a1][b1] * ( + PB_0 * PB_0 * PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_1 * (-1.0))
                            + delta[a0][b1] * delta[a0][a1] * ( + PB_0 * PB_0 * PB_1 * PA_1)
                            + 1.0 * delta[a1][b0] * ( + PB_0 * PB_1 * PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a0][a1] * ( + PB_0 * PB_1 * PB_1 * PA_1)
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1 * (-2.0))
                        )

                        + 0.25 / S1 * (
                            1.0 * ( + PB_0 * PB_0 * PA_0 * PA_0 * PA_1 * PA_1)
                            + 1.0 * ( + PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1)
                            + 1.0 * ( + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0)
                            + 1.0 * ( + PB_0 * PB_0 * PB_1 * PB_1 * PA_1 * PA_1)
                        )

                    )

                    + F4_t[2] * (

                        0.03125 / ( S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * (-2.0))
                            + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_0 * (-8.0))
                            + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * ( + PA_0 * PA_0 * (-2.0))
                            + delta[a1][b0] * delta[a1][b0] * ( + PA_0 * PA_0 * (-2.0))
                            + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * (-2.0))
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * ( + PA_0 * PA_0 * (-2.0))
                            + delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_1 * 4.0)
                            + delta[a0][b0] * delta[a1][b0] * ( + PA_0 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_1)
                            + delta[a0][b1] * delta[a1][b1] * ( + PA_0 * PA_1)
                            + delta[a0][a1] * ( + PA_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b0] * ( + PA_0 * PA_1)
                            + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_0 * PA_1)
                            + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_0 * PA_1 * (-1.0))
                            + delta[a1][b1] * delta[a0][b1] * ( + PA_0 * PA_1)
                            + delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * ( + PA_0 * PA_1 * (-1.0))
                            + 1.0 * delta[a0][b1] * delta[a1][b1] * ( + PA_0 * PA_1 * (-2.0))
                            + delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * (-1.0))
                            + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0)
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0)
                            + 1.0 * delta[a0][b0] * ( + PB_0 * PA_0 * 2.0)
                            + 1.0 * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * (-2.0))
                            + 1.0 * delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0)
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PA_0 * 2.0)
                            + delta[a0][a1] * delta[a1][b1] * ( + PB_1 * PA_0 * (-1.0))
                            + delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_0)
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PA_0 * (-2.0))
                            + delta[a0][a1] * delta[a1][b1] * ( + PB_1 * PA_0)
                            + 1.0 * delta[a0][b0] * delta[b0][b1] * ( + PB_1 * PA_0 * (-1.0))
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_0 * 6.0)
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_1 * PA_1 * (-2.0))
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_1 * PA_1 * (-8.0))
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_1 * PA_1 * (-2.0))
                            + delta[a0][b0] * delta[a0][b0] * ( + PA_1 * PA_1 * (-2.0))
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * ( + PA_1 * PA_1 * (-2.0))
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * ( + PA_1 * PA_1 * (-2.0))
                            + 1.0 * delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_1 * 2.0)
                            + 1.0 * delta[a1][b0] * ( + PB_0 * PA_1 * 2.0)
                            + 1.0 * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_1 * (-2.0))
                            + delta[a0][a1] * delta[a0][b0] * ( + PB_0 * PA_1 * (-1.0))
                            + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_1)
                            + delta[a0][a1] * delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a0][a1] * ( + PB_0 * PA_1)
                            + delta[a0][b1] * delta[a1][b0] * delta[a0][b1] * ( + PB_0 * PA_1 * 2.0)
                            + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PA_1 * 2.0)
                            + 1.0 * delta[a1][b1] * ( + PB_1 * PA_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b1] * ( + PB_1 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * ( + PB_1 * PA_1 * 6.0)
                            + delta[a0][b1] * delta[a0][a1] * ( + PB_1 * PA_1)
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PA_1 * (-2.0))
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * (-2.0))
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * (-8.0))
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0 * (-2.0))
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * (-2.0))
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                            + 1.0 * delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-1.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1)
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * 4.0)
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1 * (-1.0))
                            + 1.0 * delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-2.0))
                            + 1.0 * delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * (-2.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * (-8.0))
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1 * (-2.0))
                            + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1 * (-2.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                            + 1.0 * delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                            + 1.0 * delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                            + 1.0 * delta[a1][b1] * delta[a1][b0] * (PB_0 * PB_1)
                            + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1)
                            + delta[a0][b1] * delta[a0][b0] * (PB_0 * PB_1)
                        )

                        + 0.0625 / ( S1 * S1 ) * (
                            1.0 * ( + PA_0 * PA_0 * PA_1 * PA_1)
                            + delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * PA_1 * PA_1 * 2.0)
                            + delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * PA_0 * PA_1 * (-1.0))
                            + delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_0 * PA_0 * PA_1)
                            + delta[a1][b1] * ( + PB_1 * PA_0 * PA_0 * PA_1)
                            + 1.0 * delta[a1][b1] * ( + PB_1 * PA_0 * PA_0 * PA_1 * (-1.0))
                            + 1.0 * ( + PB_0 * PB_0 * PA_0 * PA_0)
                            + delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PB_0 * PA_0 * PA_0 * 2.0)
                            + 1.0 * ( + PB_1 * PB_1 * PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PB_1 * PA_0 * PA_0 * 2.0)
                            + delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * PA_1 * PA_1 * (-1.0))
                            + delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0 * PA_1 * PA_1)
                            + delta[a0][b1] * ( + PB_1 * PA_0 * PA_1 * PA_1)
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * PA_1 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PB_1 * PA_0 * PA_1)
                            + delta[a1][b1] * delta[a0][b0] * ( + PB_0 * PB_1 * PA_0 * PA_1 * (-1.0))
                            + delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0)
                            + 1.0 * delta[a0][b1] * ( + PB_0 * PB_0 * PB_1 * PA_0 * (-1.0))
                            + delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0)
                            + 1.0 * delta[a0][b0] * ( + PB_0 * PB_1 * PB_1 * PA_0 * (-1.0))
                            + 1.0 * ( + PB_0 * PB_0 * PA_1 * PA_1)
                            + delta[a0][b1] * delta[a0][b1] * ( + PB_0 * PB_0 * PA_1 * PA_1 * 2.0)
                            + 1.0 * ( + PB_1 * PB_1 * PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PB_1 * PA_1 * PA_1 * 2.0)
                            + 1.0 * (PB_0 * PB_0 * PB_1 * PB_1)
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1 * 2.0)
                        )

                        + 0.015625 / ( S1 * S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1]
                            + 1.0 * delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * 4.0
                            + 1.0 * delta[a1][b1] * delta[a1][b0] * delta[b0][b1]
                            + 1.0 * delta[a1][b0] * delta[a1][b0]
                            + 1.0 * delta[a1][b0] * delta[a1][b1] * delta[b0][b1]
                            + 1.0 * delta[a1][b1] * delta[a1][b1]
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * 8.0
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 19.0
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * 4.0
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 17.0
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * 2.0
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * 4.0
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1]
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0]
                            + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * 2.0
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1]
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * 7.0
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 15.0
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * 3.0
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1]
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1]
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * 3.0
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1]
                            + delta[a0][a1] * delta[a0][a1]
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0]
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1]
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * 2.0
                            + 1.0 * delta[a0][b0] * delta[a0][b0]
                            + 1.0 * delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 3.0
                            + 1.0 * delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                            + 1.0 * delta[a0][b1] * delta[a0][b1]
                            + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * 3.0
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * 7.0
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1]
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * 3.0
                        )

                    )

                    + F4_t[3] * (

                        0.0078125 / ( S1 * S1 * S1 * S1 ) * (
                            1.0 * delta[b0][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + 1.0 * delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * (1.0 * (-8.0))
                            + 1.0 * delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * (1.0 * (-2.0))
                            + 1.0 * delta[a1][b0] * delta[a1][b0] * (1.0 * (-2.0))
                            + 1.0 * delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * (1.0 * (-2.0))
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * (1.0 * (-8.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * (1.0 * (-22.0))
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * (1.0 * (-4.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (1.0 * (-8.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * (1.0 * (-18.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * (1.0 * (-4.0))
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (1.0 * (-8.0))
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * (1.0 * (-4.0))
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (1.0 * (-2.0))
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (1.0 * (-6.0))
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * (1.0 * (-14.0))
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * (1.0 * (-6.0))
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (1.0 * (-2.0))
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * (1.0 * (-2.0))
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * (1.0 * (-6.0))
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * (1.0 * (-2.0))
                            + delta[a0][a1] * delta[a0][a1] * (1.0 * (-2.0))
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (1.0 * (-2.0))
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (1.0 * (-2.0))
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * (1.0 * (-4.0))
                            + 1.0 * delta[a0][b0] * delta[a0][b0] * (1.0 * (-2.0))
                            + 1.0 * delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (1.0 * (-6.0))
                            + 1.0 * delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * (1.0 * (-2.0))
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * (1.0 * (-2.0))
                            + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * (1.0 * (-6.0))
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * (1.0 * (-6.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * (1.0 * (-2.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * (1.0 * (-6.0))
                        )

                        + 0.015625 / ( S1 * S1 * S1 ) * (
                            1.0 * ( + PA_0 * PA_0)
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[a1][b0] * ( + PA_0 * PA_0 * 2.0)
                            + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * ( + PA_0 * PA_0 * 5.0)
                            + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * ( + PA_0 * PA_0 * 2.0)
                            + delta[a1][b1] * delta[a1][b1] * ( + PA_0 * PA_0)
                            + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * ( + PA_0 * PA_0)
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * ( + PA_0 * PA_0)
                            + delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * ( + PA_0 * PA_1 * (-1.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_0 * PA_1)
                            + delta[a0][b0] * ( + PB_0 * PA_0)
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ( + PB_0 * PA_0)
                            + 1.0 * delta[a0][b0] * ( + PB_0 * PA_0 * (-1.0))
                            + 1.0 * delta[b0][b1] * delta[a0][b1] * ( + PB_0 * PA_0)
                            + 1.0 * delta[a0][b1] * delta[b0][b1] * ( + PB_0 * PA_0 * (-1.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b1] * ( + PB_0 * PA_0 * (-1.0))
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b0] * ( + PB_1 * PA_0)
                            + delta[a0][b1] * ( + PB_1 * PA_0)
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * ( + PB_1 * PA_0 * 2.0)
                            + 1.0 * delta[a0][b1] * ( + PB_1 * PA_0 * (-1.0))
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_0 * (-2.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PA_0 * (-1.0))
                            + 1.0 * ( + PA_1 * PA_1)
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * ( + PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b0] * ( + PA_1 * PA_1 * 2.0)
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ( + PA_1 * PA_1 * 5.0)
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * ( + PA_1 * PA_1 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * ( + PA_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * ( + PA_1 * PA_1)
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * ( + PA_1 * PA_1)
                            + 1.0 * delta[a1][b1] * delta[b0][b1] * ( + PB_0 * PA_1 * (-1.0))
                            + 1.0 * delta[b0][b1] * delta[a1][b1] * ( + PB_0 * PA_1)
                            + delta[a0][b0] * delta[a1][b1] * delta[a0][b1] * ( + PB_0 * PA_1)
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b1] * ( + PB_0 * PA_1 * (-1.0))
                            + 1.0 * delta[a1][b1] * ( + PB_1 * PA_1)
                            + 1.0 * delta[a1][b1] * ( + PB_1 * PA_1 * (-1.0))
                            + delta[a0][b0] * delta[a1][b1] * delta[a0][b0] * ( + PB_1 * PA_1)
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * ( + PB_1 * PA_1 * (-2.0))
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * ( + PB_1 * PA_1)
                            + delta[a0][b1] * delta[a1][b0] * delta[a0][b0] * ( + PB_1 * PA_1)
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * ( + PB_1 * PA_1 * (-1.0))
                            + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * ( + PB_1 * PA_1 * (-1.0))
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * ( + PB_1 * PA_1)
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1)
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-1.0))
                            + 1.0 * (PB_0 * PB_0)
                            + 1.0 * delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * 2.0)
                            + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * 2.0)
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0 * 2.0)
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * 4.0)
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0 * 2.0)
                            + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                            + 1.0 * delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                            + 1.0 * (PB_1 * PB_1)
                            + 1.0 * delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * 2.0)
                            + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1 * 2.0)
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1 * 2.0)
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * 4.0)
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1 * 2.0)
                            + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                            + 1.0 * delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                        )

                    )

                    + F4_t[4] * (

                        0.00390625 / ( S1 * S1 * S1 * S1 ) * (
                            1.0
                            + 1.0 * delta[b0][b1] * delta[b0][b1] * 2.0
                            + 1.0 * delta[a1][b0] * delta[a1][b0] * 2.0
                            + 1.0 * delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * 5.0
                            + 1.0 * delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * 2.0
                            + 1.0 * delta[a1][b1] * delta[a1][b1]
                            + 1.0 * delta[a1][b0] * delta[a1][b1] * delta[b0][b1]
                            + 1.0 * delta[a1][b1] * delta[a1][b1]
                            + delta[a0][a1] * delta[a0][a1] * 2.0
                            + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * 4.0
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * 2.0
                            + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 10.0
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * 3.0
                            + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * 2.0
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * 4.0
                            + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 7.0
                            + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * 2.0
                            + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * 4.0
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1]
                            + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * 2.0
                            + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * 3.0
                            + delta[a0][b0] * delta[a0][b0]
                            + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 2.0
                            + delta[a0][b0] * delta[a1][b1] * delta[a0][b0] * delta[a1][b1]
                            + delta[a0][b0] * delta[a1][b1] * delta[a1][b0] * delta[a0][b1]
                            + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * 2.0
                            + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 5.0
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * 4.0
                            + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * 2.0
                            + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                            + delta[a0][b1] * delta[a0][b1]
                            + delta[a0][b1] * delta[a1][b0] * delta[a0][b0] * delta[a1][b1]
                            + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * 2.0
                            + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * 3.0
                            + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * 2.0
                            + 1.0 * delta[a0][b0] * delta[a0][b0]
                            + 1.0 * delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 3.0
                            + 1.0 * delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                            + 1.0 * delta[a0][b1] * delta[a0][b1]
                            + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * 3.0
                            + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * 2.0
                            + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1]
                            + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * 3.0
                        )

                    )

                );

        mat_Q[ij] = sqrt(eri_ijij);
    }
}

}  // namespace gpu
