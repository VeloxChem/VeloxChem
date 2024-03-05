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
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "StringFormat.hpp"

#define TILE_DIM 16

#define MATH_CONST_PI 3.14159265358979323846

#define MATH_CONST_HALF_SQRT_PI 0.88622692545275794096

#define MATH_CONST_TWO_OVER_SQRT_PI 1.12837916709551255856

#define hipSafe(e)                                                                                                      \
    {                                                                                                                   \
        hipError_t err = (e);                                                                                           \
        if (err != hipSuccess)                                                                                          \
        {                                                                                                               \
            std::cerr << "HIP error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                    \
        }                                                                                                               \
    }

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

        const double fx = std::exp(-fa);

        for (uint32_t j = 0; j < N; j++)
        {
            values[N - j - 1] = ft[N - j - 1] * (f2a * values[N - j] + fx);
        }
    }
    else
    {
        const double fia = 1.0 / fa;

        double pf = 0.5 * fia;

        values[0] = MATH_CONST_HALF_SQRT_PI * std::sqrt(fia);

        if (pnt < 921)
        {
            const double fia2 = fia * fia;

            const double f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

            const double fx = std::exp(-fa);

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

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

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

        const auto j_cart = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto K_ij_00 = 3.0 * a_i * a_j / (a_i + a_j) - 2.0 * a_i * a_i * a_j * a_j / ((a_i + a_j) * (a_i + a_j)) * r2_ij;

        const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

        const auto K_ij_01 = -2.0 * a_i * a_i * a_j / ((a_i + a_j) * (a_i + a_j)) * rij[j_cart];

        mat_S[ij] = S_ij_01;

        mat_T[ij] = K_ij_01 * S_ij_00 + K_ij_00 * S_ij_01;
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
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        // TODO: used shared memory?
        const uint32_t d_cart_inds_0[6] = {0,0,0,1,1,2};
        const uint32_t d_cart_inds_1[6] = {0,1,2,1,2,2};

        const auto j_cart_0 = d_cart_inds_0[j % 6];
        const auto j_cart_1 = d_cart_inds_1[j % 6];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[j_cart_0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[j_cart_1];

        mat_S[ij] = (PB_0 * PB_1 + static_cast<double>(j_cart_0 == j_cart_1) * (0.5 / (a_i + a_j))) * S_ij_00;

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_T[ij] = (

                (PB_0 * PB_1 + static_cast<double>(j_cart_0 == j_cart_1) * (0.5 / (a_i + a_j))) * K_ij_00

                + 4.0 * ksi * PB_0 * PB_1

                ) * S_ij_00;
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
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        const auto i_cart = i % 3;
        const auto j_cart = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto K_ij_00 = 3.0 * a_i * a_j / (a_i + a_j) - 2.0 * a_i * a_i * a_j * a_j / ((a_i + a_j) * (a_i + a_j)) * r2_ij;

        const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;

        const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

        const auto S_ij_11 = (

                static_cast<double>(i_cart == j_cart) * (0.5 / (a_i + a_j)) +

                (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                ) * S_ij_00;

        const auto K_ij_10 = -2.0 * a_i * a_j * a_j / ((a_i + a_j) * (a_i + a_j)) * (-rij[i_cart]);

        const auto K_ij_01 = -2.0 * a_i * a_i * a_j / ((a_i + a_j) * (a_i + a_j)) * rij[j_cart];

        const auto K_ij_11 = a_i * a_j / ((a_i + a_j) * (a_i + a_j)) * static_cast<double>(i_cart == j_cart);

        mat_S[ij] = S_ij_11;

        mat_T[ij] = K_ij_11 * S_ij_00 + K_ij_10 * S_ij_01 + K_ij_01 * S_ij_10 + K_ij_00 * S_ij_11;
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
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        const auto i_cart = i % 3;

        const uint32_t d_cart_inds_0[6] = {0,0,0,1,1,2};
        const uint32_t d_cart_inds_1[6] = {0,1,2,1,2,2};

        const auto j_cart_0 = d_cart_inds_0[j % 6];
        const auto j_cart_1 = d_cart_inds_1[j % 6];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA = (a_j / (a_i + a_j)) * rij[i_cart];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[j_cart_0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[j_cart_1];

        mat_S[ij] = (PA * PB_0 * PB_1

                + (PA * static_cast<double>(j_cart_0 == j_cart_1)
                 + PB_1 * static_cast<double>(i_cart == j_cart_0)
                 + PB_0 * static_cast<double>(i_cart == j_cart_1)) * (0.5 / (a_i + a_j))

                ) * S_ij_00;

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_T[ij] = (

                (PA * PB_0 * PB_1
                 + (PA * static_cast<double>(j_cart_0 == j_cart_1)
                  + PB_0 * static_cast<double>(i_cart == j_cart_1)
                  + PB_1 * static_cast<double>(i_cart == j_cart_0)) * (0.5 / (a_i + a_j))) * K_ij_00

                + ksi * (6.0 * PA * PB_0 * PB_1
                       + (4.0 * PB_0 * static_cast<double>(i_cart == j_cart_1)
                        + 4.0 * PB_1 * static_cast<double>(i_cart == j_cart_0)
                        + 2.0 * PA * static_cast<double>(j_cart_0 == j_cart_1)) * (0.5 / (a_i + a_j)))

                ) * S_ij_00;
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
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        const uint32_t d_cart_inds_0[6] = {0,0,0,1,1,2};
        const uint32_t d_cart_inds_1[6] = {0,1,2,1,2,2};

        const auto i_cart_0 = d_cart_inds_0[i % 6];
        const auto i_cart_1 = d_cart_inds_1[i % 6];

        const auto j_cart_0 = d_cart_inds_0[j % 6];
        const auto j_cart_1 = d_cart_inds_1[j % 6];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[i_cart_0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[i_cart_1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[j_cart_0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[j_cart_1];

        //const auto K_ij_00 = 3.0 * a_i * a_j / (a_i + a_j) - 2.0 * a_i * a_i * a_j * a_j / ((a_i + a_j) * (a_i + a_j)) * r2_ij;
        //const auto K_ij_10 = -2.0 * a_i * a_j * a_j / ((a_i + a_j) * (a_i + a_j)) * (-rij[i_cart]);
        //const auto K_ij_01 = -2.0 * a_i * a_i * a_j / ((a_i + a_j) * (a_i + a_j)) * rij[j_cart];
        //const auto K_ij_11 = a_i * a_j / ((a_i + a_j) * (a_i + a_j)) * static_cast<double>(i_cart == j_cart);

        mat_S[ij] = (PA_0 * PA_1 * PB_0 * PB_1

                + (PA_1 * PA_0 * static_cast<double>(j_cart_0 == j_cart_1)
                 + PA_1 * PB_1 * static_cast<double>(i_cart_0 == j_cart_0)
                 + PA_1 * PB_0 * static_cast<double>(i_cart_0 == j_cart_1)
                 + PA_0 * PB_1 * static_cast<double>(i_cart_1 == j_cart_0)
                 + PA_0 * PB_0 * static_cast<double>(i_cart_1 == j_cart_1)
                 + PB_1 * PB_0 * static_cast<double>(i_cart_0 == i_cart_1)) * (0.5 / (a_i + a_j))

                + (static_cast<double>(i_cart_0 == i_cart_1) * static_cast<double>(j_cart_0 == j_cart_1)
                 + static_cast<double>(i_cart_0 == j_cart_1) * static_cast<double>(i_cart_1 == j_cart_0)
                 + static_cast<double>(i_cart_0 == j_cart_0) * static_cast<double>(i_cart_1 == j_cart_1)) * (0.25 / ((a_i + a_j) * (a_i + a_j)))

                ) * S_ij_00;

        const auto ksi = a_i * a_j / (a_i + a_j);

        const auto K_ij_00 = ksi * (3.0 - 2.0 * ksi * r2_ij);

        mat_T[ij] = (

                (PA_0 * PA_1 * PB_0 * PB_1

                 + (PA_1 * PA_0 * static_cast<double>(j_cart_0 == j_cart_1)
                  + PA_1 * PB_1 * static_cast<double>(i_cart_0 == j_cart_0)
                  + PA_1 * PB_0 * static_cast<double>(i_cart_0 == j_cart_1)
                  + PA_0 * PB_1 * static_cast<double>(i_cart_1 == j_cart_0)
                  + PA_0 * PB_0 * static_cast<double>(i_cart_1 == j_cart_1)
                  + PB_1 * PB_0 * static_cast<double>(i_cart_0 == i_cart_1)) * (0.5 / (a_i + a_j))

                 + (static_cast<double>(i_cart_0 == i_cart_1) * static_cast<double>(j_cart_0 == j_cart_1)
                  + static_cast<double>(i_cart_0 == j_cart_1) * static_cast<double>(i_cart_1 == j_cart_0)
                  + static_cast<double>(i_cart_0 == j_cart_0) * static_cast<double>(i_cart_1 == j_cart_1)) * (0.25 / ((a_i + a_j) * (a_i + a_j)))

                 ) * K_ij_00

                + ksi * (8.0 * PA_0 * PA_1 * PB_0 * PB_1

                       + (4.0 * PA_1 * PA_0 * static_cast<double>(j_cart_0 == j_cart_1)
                        + 6.0 * PA_1 * PB_1 * static_cast<double>(i_cart_0 == j_cart_0)
                        + 6.0 * PA_1 * PB_0 * static_cast<double>(i_cart_0 == j_cart_1)
                        + 6.0 * PA_0 * PB_1 * static_cast<double>(i_cart_1 == j_cart_0)
                        + 6.0 * PA_0 * PB_0 * static_cast<double>(i_cart_1 == j_cart_1)
                        + 4.0 * PB_1 * PB_0 * static_cast<double>(i_cart_0 == i_cart_1)) * (0.5 / (a_i + a_j))

                       + (4.0 * static_cast<double>(i_cart_0 == j_cart_1) * static_cast<double>(i_cart_1 == j_cart_0)
                        + 4.0 * static_cast<double>(i_cart_0 == j_cart_0) * static_cast<double>(i_cart_1 == j_cart_1)) * (0.25 / ((a_i + a_j) * (a_i + a_j)))
                    )

                ) * S_ij_00;
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

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * std::sqrt(a_i + a_j);

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double rPC[3] = {y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PC = rPC[0] * rPC[0] + rPC[1] * rPC[1] + rPC[2] * rPC[2];

            double F0_t[1];

            gpu::computeBoysFunction(F0_t, (a_i + a_j) * r2_PC, 0, boys_func_table, boys_func_ft);

            V_ij += q_c * S_ij_00 * F0_t[0];
        }

        mat_V[ij] = theta * V_ij;
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

        const auto j_cart = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * std::sqrt(a_i + a_j);

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double rPC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                   y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = rPC[0] * rPC[0] + rPC[1] * rPC[1] + rPC[2] * rPC[2];

            double F1_t[2];

            gpu::computeBoysFunction(F1_t, (a_i + a_j) * r2_PC, 1, boys_func_table, boys_func_ft);

            V_ij += q_c * (S_ij_01 * F1_t[0] + S_ij_00 * rPC[j_cart] * F1_t[1]);
        }

        mat_V[ij] = theta * V_ij;
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
    // each thread computes a primitive S/T matrix element

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

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

        const auto i_cart = i % 3;
        const auto j_cart = j % 3;

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;

        const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

        const auto S_ij_11 = (

                static_cast<double>(i_cart == j_cart) * (0.5 / (a_i + a_j)) +

                (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                ) * S_ij_00;

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * std::sqrt(a_i + a_j);

        double V_ij = 0.0;

        for (uint32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double rPC[3] = {x_c - (a_i * x_i + a_j * x_j) / (a_i + a_j),
                                   y_c - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   z_c - (a_i * z_i + a_j * z_j) / (a_i + a_j)};

            const auto r2_PC = rPC[0] * rPC[0] + rPC[1] * rPC[1] + rPC[2] * rPC[2];

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * r2_PC, 2, boys_func_table, boys_func_ft);

            V_ij += q_c * (S_ij_11 * F2_t[0]

                         + S_ij_10 * rPC[j_cart] * F2_t[1]

                         + S_ij_01 * rPC[i_cart] * F2_t[1]

                         + S_ij_00 * (rPC[i_cart] * rPC[j_cart] * F2_t[2] - 0.5 / (a_i + a_j) * static_cast<double>(i_cart == j_cart) * F2_t[1]));
        }

        mat_V[ij] = theta * V_ij;
    }
}

}  // namespace gpu
