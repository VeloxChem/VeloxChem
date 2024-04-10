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

#include "AngularMomentumIntegrals.hpp"
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
computeAngularMomentumSS(double*         mat_mu_X,
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

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            /*
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            mat_mu[m][ij] = S_ij_00 * L_ij;
            */

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2]
                        )

                    )

                    + (-1) * PC[m1] * (

                        -1.0 * 2.0 * a_j * (
                            PB_m2
                        )

                    )

                    + (-1) * PC[m2] * (

                        2.0 * a_j * (
                            PB_m1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeAngularMomentumSP(double*         mat_mu_X,
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

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const double one_CA[3][3] = {{0.0,   -CA[2], CA[1]},
                                     {CA[2],  0.0,  -CA[0]},
                                     {-CA[1], CA[0],  0.0}};

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PB_0)
                            + delta[b0][m1] * (PB_m2)
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PB_0)
                            + delta[b0][m2] * (PB_m1)
                        )

                    )

                    + (-1) * PC[m1] * (

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b0][m2]
                        )

                        -1.0 * 2.0 * a_j * (
                            PB_0 * PB_m2
                        )

                        + (
                            delta[b0][m2]
                        )

                    )

                    + (-1) * PC[m2] * (

                        -1 * (
                            delta[b0][m1]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b0][m1]
                        )

                        + 2.0 * a_j * (
                            PB_0 * PB_m1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeAngularMomentumSD(double*         mat_mu_X,
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

        const double PC[3] = {origin_X - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                              origin_Y - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                              origin_Z - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const double one_CA[3][3] = {{0.0,   -CA[2], CA[1]},
                                     {CA[2],  0.0,  -CA[0]},
                                     {-CA[1], CA[0],  0.0}};

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2]
                            + delta[b0][m1] * delta[b1][m2]
                            + delta[b0][m2] * delta[b1][m1]
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PB_0 * PB_1)
                            + delta[b1][m1] * (PB_0 * PB_m2)
                            + delta[b0][m1] * (PB_1 * PB_m2)
                        )

                        -0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2]
                            + delta[b0][m2] * delta[b1][m1]
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2]
                            + delta[b0][m1] * delta[b1][m2]
                            + delta[b0][m2] * delta[b1][m1]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PB_0 * PB_1)
                            + delta[b1][m2] * (PB_0 * PB_m1)
                            + delta[b0][m2] * (PB_1 * PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2]
                            + delta[b0][m2] * delta[b1][m1]
                        )

                    )

                    + (-1) * PC[m1] * (

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m2] * (PB_0)
                            + delta[b0][m2] * (PB_1)
                            + delta[b0][b1] * (PB_m2)
                        )

                        -1.0 * 2.0 * a_j * (
                            PB_0 * PB_1 * PB_m2
                        )

                        + (
                            delta[b1][m2] * (PB_0)
                            + delta[b0][m2] * (PB_1)
                        )

                    )

                    + (-1) * PC[m2] * (

                        -1 * (
                            delta[b1][m1] * (PB_0)
                            + delta[b0][m1] * (PB_1)
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m1] * (PB_0)
                            + delta[b0][m1] * (PB_1)
                            + delta[b0][b1] * (PB_m1)
                        )

                        + 2.0 * a_j * (
                            PB_0 * PB_1 * PB_m1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeAngularMomentumPP(double*         mat_mu_X,
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
    __shared__ double   one_one[3][3][3];

    // each thread computes a primitive matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        // x cross x/y/z
        one_one[0][0][0] =  0.0; one_one[0][0][1] =  0.0; one_one[0][0][2] =  0.0;
        one_one[0][1][0] =  0.0; one_one[0][1][1] =  0.0; one_one[0][1][2] =  1.0;
        one_one[0][2][0] =  0.0; one_one[0][2][1] = -1.0; one_one[0][2][2] =  0.0;

        // y cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  0.0; one_one[1][0][2] = -1.0;
        one_one[1][1][0] =  0.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  1.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;

        // z cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  1.0; one_one[1][0][2] =  0.0;
        one_one[1][1][0] = -1.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  0.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;
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

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const double one_CA[3][3] = {{0.0,   -CA[2], CA[1]},
                                     {CA[2],  0.0,  -CA[0]},
                                     {-CA[1], CA[0],  0.0}};

        const double one_CB[3][3] = {{0.0,   -CB[2], CB[1]},
                                     {CB[2],  0.0,  -CB[0]},
                                     {-CB[1], CB[0],  0.0}};

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PA_m1 = (a_j / (a_i + a_j)) * rij[m1];
            const auto PA_m2 = (a_j / (a_i + a_j)) * rij[m2];

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][b0] * delta[m1][m2]
                            + delta[a0][m1] * delta[b0][m2]
                            + delta[a0][m2] * delta[b0][m1]
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PB_0)
                            + delta[b0][m1] * (PA_0 * PB_m2)
                            + delta[a0][m1] * (PB_0 * PB_m2)
                        )

                        -0.5 / (a_i + a_j) * (
                            delta[a0][m2] * delta[b0][m1]
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][b0] * delta[m1][m2]
                            + delta[a0][m1] * delta[b0][m2]
                            + delta[a0][m2] * delta[b0][m1]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PB_0)
                            + delta[b0][m2] * (PA_0 * PB_m1)
                            + delta[a0][m2] * (PB_0 * PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[a0][m1] * delta[b0][m2]
                        )

                    )

                    + (-1) * PC[m1] * (

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b0][m2] * (PA_0)
                            + delta[a0][m2] * (PB_0)
                            + delta[a0][b0] * (PB_m2)
                        )

                        -1.0 * 2.0 * a_j * (
                            PA_0 * PB_0 * PB_m2
                        )

                        + (
                            delta[b0][m2] * (PA_0)
                        )

                    )

                    + (-1) * PC[m2] * (

                        -1 * (
                            delta[b0][m1] * (PA_0)
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b0][m1] * (PA_0)
                            + delta[a0][m1] * (PB_0)
                            + delta[a0][b0] * (PB_m1)
                        )

                        + 2.0 * a_j * (
                            PA_0 * PB_0 * PB_m1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeAngularMomentumPD(double*         mat_mu_X,
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
    __shared__ double   one_one[3][3][3];

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

        // x cross x/y/z
        one_one[0][0][0] =  0.0; one_one[0][0][1] =  0.0; one_one[0][0][2] =  0.0;
        one_one[0][1][0] =  0.0; one_one[0][1][1] =  0.0; one_one[0][1][2] =  1.0;
        one_one[0][2][0] =  0.0; one_one[0][2][1] = -1.0; one_one[0][2][2] =  0.0;

        // y cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  0.0; one_one[1][0][2] = -1.0;
        one_one[1][1][0] =  0.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  1.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;

        // z cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  1.0; one_one[1][0][2] =  0.0;
        one_one[1][1][0] = -1.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  0.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;
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

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const double one_CA[3][3] = {{0.0,   -CA[2], CA[1]},
                                     {CA[2],  0.0,  -CA[0]},
                                     {-CA[1], CA[0],  0.0}};

        const double one_CB[3][3] = {{0.0,   -CB[2], CB[1]},
                                     {CB[2],  0.0,  -CB[0]},
                                     {-CB[1], CB[0],  0.0}};

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PA_m1 = (a_j / (a_i + a_j)) * rij[m1];
            const auto PA_m2 = (a_j / (a_i + a_j)) * rij[m2];

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2] * (PA_0)
                            + delta[b0][m1] * delta[b1][m2] * (PA_0)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0)
                            + delta[a0][b1] * delta[m1][m2] * (PB_0)
                            + delta[a0][m1] * delta[b1][m2] * (PB_0)
                            + delta[a0][m2] * delta[b1][m1] * (PB_0)
                            + delta[a0][b0] * delta[m1][m2] * (PB_1)
                            + delta[a0][m1] * delta[b0][m2] * (PB_1)
                            + delta[a0][m2] * delta[b0][m1] * (PB_1)
                            + delta[a0][b0] * delta[b1][m1] * (PB_m2)
                            + delta[a0][b1] * delta[b0][m1] * (PB_m2)
                            + delta[a0][m1] * delta[b0][b1] * (PB_m2)
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PB_0 * PB_1)
                            + delta[b1][m1] * (PA_0 * PB_0 * PB_m2)
                            + delta[b0][m1] * (PA_0 * PB_1 * PB_m2)
                            + delta[a0][m1] * (PB_0 * PB_1 * PB_m2)
                        )

                        -0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2] * (PA_0)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0)
                            + delta[a0][m2] * delta[b1][m1] * (PB_0)
                            + delta[a0][m2] * delta[b0][m1] * (PB_1)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2] * (PA_0)
                            + delta[b0][m1] * delta[b1][m2] * (PA_0)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0)
                            + delta[a0][b1] * delta[m1][m2] * (PB_0)
                            + delta[a0][m1] * delta[b1][m2] * (PB_0)
                            + delta[a0][m2] * delta[b1][m1] * (PB_0)
                            + delta[a0][b0] * delta[m1][m2] * (PB_1)
                            + delta[a0][m1] * delta[b0][m2] * (PB_1)
                            + delta[a0][m2] * delta[b0][m1] * (PB_1)
                            + delta[a0][b0] * delta[b1][m2] * (PB_m1)
                            + delta[a0][b1] * delta[b0][m2] * (PB_m1)
                            + delta[a0][m2] * delta[b0][b1] * (PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PB_0 * PB_1)
                            + delta[b1][m2] * (PA_0 * PB_0 * PB_m1)
                            + delta[b0][m2] * (PA_0 * PB_1 * PB_m1)
                            + delta[a0][m2] * (PB_0 * PB_1 * PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2] * (PA_0)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0)
                            + delta[a0][m1] * delta[b1][m2] * (PB_0)
                            + delta[a0][m1] * delta[b0][m2] * (PB_1)
                        )

                    )

                    + (-1) * PC[m1] * (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][b0] * delta[b1][m2]
                            + delta[a0][b1] * delta[b0][m2]
                            + delta[a0][m2] * delta[b0][b1]
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m2] * (PA_0 * PB_0)
                            + delta[b0][m2] * (PA_0 * PB_1)
                            + delta[b0][b1] * (PA_0 * PB_m2)
                            + delta[a0][m2] * (PB_0 * PB_1)
                            + delta[a0][b1] * (PB_0 * PB_m2)
                            + delta[a0][b0] * (PB_1 * PB_m2)
                        )

                        -1.0 * 2.0 * a_j * (
                            PA_0 * PB_0 * PB_1 * PB_m2
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[a0][b0] * delta[b1][m2]
                            + delta[a0][b1] * delta[b0][m2]
                        )

                        + (
                            delta[b1][m2] * (PA_0 * PB_0)
                            + delta[b0][m2] * (PA_0 * PB_1)
                        )

                    )

                    + (-1) * PC[m2] * (

                        -0.5 / (a_i + a_j) * (
                            delta[a0][b0] * delta[b1][m1]
                            + delta[a0][b1] * delta[b0][m1]
                        )

                        -1 * (
                            delta[b1][m1] * (PA_0 * PB_0)
                            + delta[b0][m1] * (PA_0 * PB_1)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][b0] * delta[b1][m1]
                            + delta[a0][b1] * delta[b0][m1]
                            + delta[a0][m1] * delta[b0][b1]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m1] * (PA_0 * PB_0)
                            + delta[b0][m1] * (PA_0 * PB_1)
                            + delta[b0][b1] * (PA_0 * PB_m1)
                            + delta[a0][m1] * (PB_0 * PB_1)
                            + delta[a0][b1] * (PB_0 * PB_m1)
                            + delta[a0][b0] * (PB_1 * PB_m1)
                        )

                        + 2.0 * a_j * (
                            PA_0 * PB_0 * PB_1 * PB_m1
                        )

                    )

                    );
        }
    }
}

__global__ void
computeAngularMomentumDD(double*         mat_mu_X,
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
    __shared__ double   one_one[3][3][3];

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

        // x cross x/y/z
        one_one[0][0][0] =  0.0; one_one[0][0][1] =  0.0; one_one[0][0][2] =  0.0;
        one_one[0][1][0] =  0.0; one_one[0][1][1] =  0.0; one_one[0][1][2] =  1.0;
        one_one[0][2][0] =  0.0; one_one[0][2][1] = -1.0; one_one[0][2][2] =  0.0;

        // y cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  0.0; one_one[1][0][2] = -1.0;
        one_one[1][1][0] =  0.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  1.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;

        // z cross x/y/z
        one_one[1][0][0] =  0.0; one_one[1][0][1] =  1.0; one_one[1][0][2] =  0.0;
        one_one[1][1][0] = -1.0; one_one[1][1][1] =  0.0; one_one[1][1][2] =  0.0;
        one_one[1][2][0] =  0.0; one_one[1][2][1] =  0.0; one_one[1][2][2] =  0.0;
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

        const double CA[3] = {x_i - origin_X, y_i - origin_Y, z_i - origin_Z};
        const double CB[3] = {x_j - origin_X, y_j - origin_Y, z_j - origin_Z};

        const double CA_CB[3] = {
            CA[1] * CB[2] - CA[2] * CB[1],
            CA[2] * CB[0] - CA[0] * CB[2],
            CA[0] * CB[1] - CA[1] * CB[0],
        };

        const double one_CA[3][3] = {{0.0,   -CA[2], CA[1]},
                                     {CA[2],  0.0,  -CA[0]},
                                     {-CA[1], CA[0],  0.0}};

        const double one_CB[3][3] = {{0.0,   -CB[2], CB[1]},
                                     {CB[2],  0.0,  -CB[0]},
                                     {-CB[1], CB[0],  0.0}};

        const auto ksi = a_i * a_j / (a_i + a_j);

        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[3] = {mat_mu_X, mat_mu_Y, mat_mu_Z};

        for (uint32_t m = 0; m < 3; m++)
        {
            const auto L_ij = 2.0 * ksi * CA_CB[m];

            const auto m1 = (m + 1) % 3;
            const auto m2 = (m + 2) % 3;

            const auto PA_m1 = (a_j / (a_i + a_j)) * rij[m1];
            const auto PA_m2 = (a_j / (a_i + a_j)) * rij[m2];

            const auto PB_m1 = (-a_i / (a_i + a_j)) * rij[m1];
            const auto PB_m2 = (-a_i / (a_i + a_j)) * rij[m2];

            mat_mu[m][ij] = S_ij_00 * (

                    (

                        -0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][a1] * delta[b0][b1] * delta[m1][m2]
                            + delta[a0][a1] * delta[b0][m1] * delta[b1][m2]
                            + delta[a0][a1] * delta[b0][m2] * delta[b1][m1]
                            + delta[a0][b0] * delta[a1][b1] * delta[m1][m2]
                            + delta[a0][b0] * delta[a1][m1] * delta[b1][m2]
                            + delta[a0][b0] * delta[a1][m2] * delta[b1][m1]
                            + delta[a0][b1] * delta[a1][b0] * delta[m1][m2]
                            + delta[a0][b1] * delta[a1][m1] * delta[b0][m2]
                            + delta[a0][b1] * delta[a1][m2] * delta[b0][m1]
                            + delta[a0][m1] * delta[a1][b0] * delta[b1][m2]
                            + delta[a0][m1] * delta[a1][b1] * delta[b0][m2]
                            + delta[a0][m1] * delta[a1][m2] * delta[b0][b1]
                            + delta[a0][m2] * delta[a1][b0] * delta[b1][m1]
                            + delta[a0][m2] * delta[a1][b1] * delta[b0][m1]
                            + delta[a0][m2] * delta[a1][m1] * delta[b0][b1]
                        )

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2] * (PA_0 * PA_1)
                            + delta[b0][m1] * delta[b1][m2] * (PA_0 * PA_1)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0 * PA_1)
                            + delta[a1][b1] * delta[m1][m2] * (PA_0 * PB_0)
                            + delta[a1][m1] * delta[b1][m2] * (PA_0 * PB_0)
                            + delta[a1][m2] * delta[b1][m1] * (PA_0 * PB_0)
                            + delta[a1][b0] * delta[m1][m2] * (PA_0 * PB_1)
                            + delta[a1][m1] * delta[b0][m2] * (PA_0 * PB_1)
                            + delta[a1][m2] * delta[b0][m1] * (PA_0 * PB_1)
                            + delta[a1][b0] * delta[b1][m1] * (PA_0 * PB_m2)
                            + delta[a1][b1] * delta[b0][m1] * (PA_0 * PB_m2)
                            + delta[a1][m1] * delta[b0][b1] * (PA_0 * PB_m2)
                            + delta[a0][b1] * delta[m1][m2] * (PA_1 * PB_0)
                            + delta[a0][m1] * delta[b1][m2] * (PA_1 * PB_0)
                            + delta[a0][m2] * delta[b1][m1] * (PA_1 * PB_0)
                            + delta[a0][b0] * delta[m1][m2] * (PA_1 * PB_1)
                            + delta[a0][m1] * delta[b0][m2] * (PA_1 * PB_1)
                            + delta[a0][m2] * delta[b0][m1] * (PA_1 * PB_1)
                            + delta[a0][b0] * delta[b1][m1] * (PA_1 * PB_m2)
                            + delta[a0][b1] * delta[b0][m1] * (PA_1 * PB_m2)
                            + delta[a0][m1] * delta[b0][b1] * (PA_1 * PB_m2)
                            + delta[a0][a1] * delta[m1][m2] * (PB_0 * PB_1)
                            + delta[a0][m1] * delta[a1][m2] * (PB_0 * PB_1)
                            + delta[a0][m2] * delta[a1][m1] * (PB_0 * PB_1)
                            + delta[a0][a1] * delta[b1][m1] * (PB_0 * PB_m2)
                            + delta[a0][b1] * delta[a1][m1] * (PB_0 * PB_m2)
                            + delta[a0][m1] * delta[a1][b1] * (PB_0 * PB_m2)
                            + delta[a0][a1] * delta[b0][m1] * (PB_1 * PB_m2)
                            + delta[a0][b0] * delta[a1][m1] * (PB_1 * PB_m2)
                            + delta[a0][m1] * delta[a1][b0] * (PB_1 * PB_m2)
                        )

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[b0][m1] * delta[b1][m2]
                            + delta[a0][a1] * delta[b0][m2] * delta[b1][m1]
                            + delta[a0][b0] * delta[a1][m2] * delta[b1][m1]
                            + delta[a0][b1] * delta[a1][m2] * delta[b0][m1]
                            + delta[a0][m2] * delta[a1][b0] * delta[b1][m1]
                            + delta[a0][m2] * delta[a1][b1] * delta[b0][m1]
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PA_1 * PB_0 * PB_1)
                            + delta[b1][m1] * (PA_0 * PA_1 * PB_0 * PB_m2)
                            + delta[b0][m1] * (PA_0 * PA_1 * PB_1 * PB_m2)
                            + delta[a1][m1] * (PA_0 * PB_0 * PB_1 * PB_m2)
                            + delta[a0][m1] * (PA_1 * PB_0 * PB_1 * PB_m2)
                        )

                        -0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2] * (PA_0 * PA_1)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0 * PA_1)
                            + delta[a1][m2] * delta[b1][m1] * (PA_0 * PB_0)
                            + delta[a1][m2] * delta[b0][m1] * (PA_0 * PB_1)
                            + delta[a0][m2] * delta[b1][m1] * (PA_1 * PB_0)
                            + delta[a0][m2] * delta[b0][m1] * (PA_1 * PB_1)
                        )

                        + 0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a0][a1] * delta[b0][b1] * delta[m1][m2]
                            + delta[a0][a1] * delta[b0][m1] * delta[b1][m2]
                            + delta[a0][a1] * delta[b0][m2] * delta[b1][m1]
                            + delta[a0][b0] * delta[a1][b1] * delta[m1][m2]
                            + delta[a0][b0] * delta[a1][m1] * delta[b1][m2]
                            + delta[a0][b0] * delta[a1][m2] * delta[b1][m1]
                            + delta[a0][b1] * delta[a1][b0] * delta[m1][m2]
                            + delta[a0][b1] * delta[a1][m1] * delta[b0][m2]
                            + delta[a0][b1] * delta[a1][m2] * delta[b0][m1]
                            + delta[a0][m1] * delta[a1][b0] * delta[b1][m2]
                            + delta[a0][m1] * delta[a1][b1] * delta[b0][m2]
                            + delta[a0][m1] * delta[a1][m2] * delta[b0][b1]
                            + delta[a0][m2] * delta[a1][b0] * delta[b1][m1]
                            + delta[a0][m2] * delta[a1][b1] * delta[b0][m1]
                            + delta[a0][m2] * delta[a1][m1] * delta[b0][b1]
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[b0][b1] * delta[m1][m2] * (PA_0 * PA_1)
                            + delta[b0][m1] * delta[b1][m2] * (PA_0 * PA_1)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0 * PA_1)
                            + delta[a1][b1] * delta[m1][m2] * (PA_0 * PB_0)
                            + delta[a1][m1] * delta[b1][m2] * (PA_0 * PB_0)
                            + delta[a1][m2] * delta[b1][m1] * (PA_0 * PB_0)
                            + delta[a1][b0] * delta[m1][m2] * (PA_0 * PB_1)
                            + delta[a1][m1] * delta[b0][m2] * (PA_0 * PB_1)
                            + delta[a1][m2] * delta[b0][m1] * (PA_0 * PB_1)
                            + delta[a1][b0] * delta[b1][m2] * (PA_0 * PB_m1)
                            + delta[a1][b1] * delta[b0][m2] * (PA_0 * PB_m1)
                            + delta[a1][m2] * delta[b0][b1] * (PA_0 * PB_m1)
                            + delta[a0][b1] * delta[m1][m2] * (PA_1 * PB_0)
                            + delta[a0][m1] * delta[b1][m2] * (PA_1 * PB_0)
                            + delta[a0][m2] * delta[b1][m1] * (PA_1 * PB_0)
                            + delta[a0][b0] * delta[m1][m2] * (PA_1 * PB_1)
                            + delta[a0][m1] * delta[b0][m2] * (PA_1 * PB_1)
                            + delta[a0][m2] * delta[b0][m1] * (PA_1 * PB_1)
                            + delta[a0][b0] * delta[b1][m2] * (PA_1 * PB_m1)
                            + delta[a0][b1] * delta[b0][m2] * (PA_1 * PB_m1)
                            + delta[a0][m2] * delta[b0][b1] * (PA_1 * PB_m1)
                            + delta[a0][a1] * delta[m1][m2] * (PB_0 * PB_1)
                            + delta[a0][m1] * delta[a1][m2] * (PB_0 * PB_1)
                            + delta[a0][m2] * delta[a1][m1] * (PB_0 * PB_1)
                            + delta[a0][a1] * delta[b1][m2] * (PB_0 * PB_m1)
                            + delta[a0][b1] * delta[a1][m2] * (PB_0 * PB_m1)
                            + delta[a0][m2] * delta[a1][b1] * (PB_0 * PB_m1)
                            + delta[a0][a1] * delta[b0][m2] * (PB_1 * PB_m1)
                            + delta[a0][b0] * delta[a1][m2] * (PB_1 * PB_m1)
                            + delta[a0][m2] * delta[a1][b0] * (PB_1 * PB_m1)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            delta[a0][a1] * delta[b0][m1] * delta[b1][m2]
                            + delta[a0][a1] * delta[b0][m2] * delta[b1][m1]
                            + delta[a0][b0] * delta[a1][m1] * delta[b1][m2]
                            + delta[a0][b1] * delta[a1][m1] * delta[b0][m2]
                            + delta[a0][m1] * delta[a1][b0] * delta[b1][m2]
                            + delta[a0][m1] * delta[a1][b1] * delta[b0][m2]
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[m1][m2] * (PA_0 * PA_1 * PB_0 * PB_1)
                            + delta[b1][m2] * (PA_0 * PA_1 * PB_0 * PB_m1)
                            + delta[b0][m2] * (PA_0 * PA_1 * PB_1 * PB_m1)
                            + delta[a1][m2] * (PA_0 * PB_0 * PB_1 * PB_m1)
                            + delta[a0][m2] * (PA_1 * PB_0 * PB_1 * PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][m1] * delta[b1][m2] * (PA_0 * PA_1)
                            + delta[b0][m2] * delta[b1][m1] * (PA_0 * PA_1)
                            + delta[a1][m1] * delta[b1][m2] * (PA_0 * PB_0)
                            + delta[a1][m1] * delta[b0][m2] * (PA_0 * PB_1)
                            + delta[a0][m1] * delta[b1][m2] * (PA_1 * PB_0)
                            + delta[a0][m1] * delta[b0][m2] * (PA_1 * PB_1)
                        )

                    )

                    + (-1) * PC[m1] * (

                        -0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a1][b0] * delta[b1][m2] * (PA_0)
                            + delta[a1][b1] * delta[b0][m2] * (PA_0)
                            + delta[a1][m2] * delta[b0][b1] * (PA_0)
                            + delta[a0][b0] * delta[b1][m2] * (PA_1)
                            + delta[a0][b1] * delta[b0][m2] * (PA_1)
                            + delta[a0][m2] * delta[b0][b1] * (PA_1)
                            + delta[a0][a1] * delta[b1][m2] * (PB_0)
                            + delta[a0][b1] * delta[a1][m2] * (PB_0)
                            + delta[a0][m2] * delta[a1][b1] * (PB_0)
                            + delta[a0][a1] * delta[b0][m2] * (PB_1)
                            + delta[a0][b0] * delta[a1][m2] * (PB_1)
                            + delta[a0][m2] * delta[a1][b0] * (PB_1)
                            + delta[a0][a1] * delta[b0][b1] * (PB_m2)
                            + delta[a0][b0] * delta[a1][b1] * (PB_m2)
                            + delta[a0][b1] * delta[a1][b0] * (PB_m2)
                        )

                        -0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m2] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m2] * (PA_0 * PA_1 * PB_1)
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_m2)
                            + delta[a1][m2] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_m2)
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_m2)
                            + delta[a0][m2] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_m2)
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_m2)
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_m2)
                        )

                        -1.0 * 2.0 * a_j * (
                            PA_0 * PA_1 * PB_0 * PB_1 * PB_m2
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[a1][b0] * delta[b1][m2] * (PA_0)
                            + delta[a1][b1] * delta[b0][m2] * (PA_0)
                            + delta[a0][b0] * delta[b1][m2] * (PA_1)
                            + delta[a0][b1] * delta[b0][m2] * (PA_1)
                            + delta[a0][a1] * delta[b1][m2] * (PB_0)
                            + delta[a0][a1] * delta[b0][m2] * (PB_1)
                        )

                        + (
                            delta[b1][m2] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m2] * (PA_0 * PA_1 * PB_1)
                        )

                    )

                    + (-1) * PC[m2] * (

                        -0.5 / (a_i + a_j) * (
                            delta[a1][b0] * delta[b1][m1] * (PA_0)
                            + delta[a1][b1] * delta[b0][m1] * (PA_0)
                            + delta[a0][b0] * delta[b1][m1] * (PA_1)
                            + delta[a0][b1] * delta[b0][m1] * (PA_1)
                            + delta[a0][a1] * delta[b1][m1] * (PB_0)
                            + delta[a0][a1] * delta[b0][m1] * (PB_1)
                        )

                        -1 * (
                            delta[b1][m1] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m1] * (PA_0 * PA_1 * PB_1)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * 2.0 * a_j * (
                            delta[a1][b0] * delta[b1][m1] * (PA_0)
                            + delta[a1][b1] * delta[b0][m1] * (PA_0)
                            + delta[a1][m1] * delta[b0][b1] * (PA_0)
                            + delta[a0][b0] * delta[b1][m1] * (PA_1)
                            + delta[a0][b1] * delta[b0][m1] * (PA_1)
                            + delta[a0][m1] * delta[b0][b1] * (PA_1)
                            + delta[a0][a1] * delta[b1][m1] * (PB_0)
                            + delta[a0][b1] * delta[a1][m1] * (PB_0)
                            + delta[a0][m1] * delta[a1][b1] * (PB_0)
                            + delta[a0][a1] * delta[b0][m1] * (PB_1)
                            + delta[a0][b0] * delta[a1][m1] * (PB_1)
                            + delta[a0][m1] * delta[a1][b0] * (PB_1)
                            + delta[a0][a1] * delta[b0][b1] * (PB_m1)
                            + delta[a0][b0] * delta[a1][b1] * (PB_m1)
                            + delta[a0][b1] * delta[a1][b0] * (PB_m1)
                        )

                        + 0.5 / (a_i + a_j) * 2.0 * a_j * (
                            delta[b1][m1] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m1] * (PA_0 * PA_1 * PB_1)
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_m1)
                            + delta[a1][m1] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_m1)
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_m1)
                            + delta[a0][m1] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_m1)
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_m1)
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_m1)
                        )

                        + 2.0 * a_j * (
                            PA_0 * PA_1 * PB_0 * PB_1 * PB_m1
                        )

                    )

                    );
        }
    }
}

}  // namespace gpu
