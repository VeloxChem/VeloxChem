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

#include "ElectricDipoleIntegrals.hpp"
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

__global__ void
computeElectricDipoleSS(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
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

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (-1) * PC[m];
        }
    }
}

__global__ void
computeElectricDipoleSP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const uint32_t* sp_first_inds_local,
                        const uint32_t* sp_second_inds_local,
                        const uint32_t  sp_prim_pair_count_local)
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

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (

                    0.5 / (a_i + a_j) * delta[m][b0]

                    + (-1) * PC[m] * PB_0

                    );
        }
    }
}

__global__ void
computeElectricDipoleSD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
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

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (

                    (

                        0.5 / (a_i + a_j) * (
                            delta[m][b1] * (PB_0)
                            + delta[m][b0] * (PB_1)
                        )

                    )

                    + (-1) * PC[m] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][b1]
                        )

                        + (
                            PB_0 * PB_1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeElectricDipolePP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
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

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (

                    (

                        0.5 / (a_i + a_j) * (
                            delta[m][a0] * (PB_0)
                            + delta[m][b0] * (PA_0)
                        )

                    )

                    + (-1) * PC[m] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][a0]
                        )

                        + (
                            PB_0 * PA_0
                        )

                    )

                    );
        }
    }
}

__global__ void
computeElectricDipolePD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
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

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI / (a_i + a_j), 1.5) * exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[b0][a0] * delta[m][b1]
                            + delta[b0][b1] * delta[m][a0]
                            + delta[b1][a0] * delta[m][b0]
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[m][a0] * (PB_0 * PB_1)
                            + delta[m][b1] * (PB_0 * PA_0)
                            + delta[m][b0] * (PB_1 * PA_0)
                        )

                    )

                    + (-1) * PC[m] * (

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
}

__global__ void
computeElectricDipoleDD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
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

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            mat_mu[m][ij] = S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[m][b1] * (PB_0)
                            + delta[b1][a0] * delta[m][a1] * (PB_0)
                            + delta[b1][a1] * delta[m][a0] * (PB_0)
                            + delta[a0][a1] * delta[m][b0] * (PB_1)
                            + delta[b0][a0] * delta[m][a1] * (PB_1)
                            + delta[b0][a1] * delta[m][a0] * (PB_1)
                            + delta[b0][a1] * delta[m][b1] * (PA_0)
                            + delta[b0][b1] * delta[m][a1] * (PA_0)
                            + delta[b1][a1] * delta[m][b0] * (PA_0)
                            + delta[b0][a0] * delta[m][b1] * (PA_1)
                            + delta[b0][b1] * delta[m][a0] * (PA_1)
                            + delta[b1][a0] * delta[m][b0] * (PA_1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[m][a1] * (PB_0 * PB_1 * PA_0)
                            + delta[m][a0] * (PB_0 * PB_1 * PA_1)
                            + delta[m][b1] * (PB_0 * PA_0 * PA_1)
                            + delta[m][b0] * (PB_1 * PA_0 * PA_1)
                        )

                    )

                    + (-1) * PC[m] * (

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
}

}  // namespace gpu
