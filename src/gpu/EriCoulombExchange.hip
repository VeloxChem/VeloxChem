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
#include "EriCoulombExchange.hpp"
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

#define hipSafe(e)                                                                                                        \
    {                                                                                                                     \
        hipError_t err = (e);                                                                                             \
        if (err != hipSuccess)                                                                                            \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl;  \
            std::exit(EXIT_FAILURE);                                                                                      \
        }                                                                                                                 \
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
computeCoulombFockSSSS(double*         mat_J,
          const double*   s_prim_info,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D,
          const double*   mat_Q,
          const uint32_t* first_inds,
          const uint32_t* second_inds,
          const uint32_t  ss_prim_pair_count,
          const double*   boys_func_table,
          const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < ss_prim_pair_count) && (ij < ss_prim_pair_count) && (std::fabs(mat_Q[ij] * mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = first_inds[ij];
            const auto j = second_inds[ij];

            const auto k = first_inds[kl];
            const auto l = second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = s_prim_info[l + s_prim_count * 0];
            const auto c_l = s_prim_info[l + s_prim_count * 1];
            const auto x_l = s_prim_info[l + s_prim_count * 2];
            const auto y_l = s_prim_info[l + s_prim_count * 3];
            const auto z_l = s_prim_info[l + s_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            const auto x_PQ = (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j);
            const auto y_PQ = (a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j);
            const auto z_PQ = (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j);

            const auto r2_PQ = x_PQ * x_PQ + y_PQ * y_PQ + z_PQ * z_PQ;

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F0_t[1];

            gpu::computeBoysFunction(F0_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 0, boys_func_table, boys_func_ft);

            const auto S_ij = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = Lambda * S_ij * S_kl * F0_t[0] * mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockSSSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   ss_mat_Q,
                       const double*   sp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < sp_prim_pair_count) && (ij < ss_prim_pair_count) && (std::fabs(ss_mat_Q[ij] * sp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            const auto k = sp_first_inds[kl];
            const auto l = sp_second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F1_t[2];

            gpu::computeBoysFunction(F1_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 1, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_ij_00 * S_kl_01 * F1_t[0] +

                    S_ij_00 * S_kl_00 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[l_cart] * F1_t[1]

                    );

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * sp_mat_D[kl] * 2.0;
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockSSPP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   ss_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < pp_prim_pair_count) && (ij < ss_prim_pair_count) && (std::fabs(ss_mat_Q[ij] * pp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            const auto k = pp_first_inds[kl];
            const auto l = pp_second_inds[kl];

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

            const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];
            const auto c_k = p_prim_info[k / 3 + p_prim_count * 1];
            const auto x_k = p_prim_info[k / 3 + p_prim_count * 2];
            const auto y_k = p_prim_info[k / 3 + p_prim_count * 3];
            const auto z_k = p_prim_info[k / 3 + p_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto k_cart = k % 3;
            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 2, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};

            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_kl_10 = (a_l / (a_k + a_l)) * rkl[k_cart] * S_kl_00;

            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const auto S_kl_11 = (

                    (k_cart == l_cart ? 1.0 : 0.0) * (0.5 / (a_k + a_l)) +

                    (a_k * a_l / ((a_k + a_l) * (a_k + a_l))) * (-rkl[k_cart]) * rkl[l_cart]

                    ) * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_ij_00 * S_kl_11 * F2_t[0] +

                    S_ij_00 * S_kl_10 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[l_cart] * F2_t[1] +

                    S_ij_00 * S_kl_01 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[k_cart] * F2_t[1] +

                    S_ij_00 * S_kl_00 * ((a_i + a_j) / (a_i + a_j + a_k + a_l)) * (

                        ((a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[k_cart] * rPQ[l_cart] * F2_t[2] -

                        (k_cart == l_cart ? 1.0 : 0.0) * (0.5 / (a_k + a_l)) * F2_t[1]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * pp_mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockSPSS(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   ss_mat_D,
                       const double*   ss_mat_Q,
                       const double*   sp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < ss_prim_pair_count) && (ij < sp_prim_pair_count) && (std::fabs(sp_mat_Q[ij] * ss_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto k = ss_first_inds[kl];
            const auto l = ss_second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = s_prim_info[l + s_prim_count * 0];
            const auto c_l = s_prim_info[l + s_prim_count * 1];
            const auto x_l = s_prim_info[l + s_prim_count * 2];
            const auto y_l = s_prim_info[l + s_prim_count * 3];
            const auto z_l = s_prim_info[l + s_prim_count * 4];

            const auto j_cart = j % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F1_t[2];

            gpu::computeBoysFunction(F1_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 1, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};

            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

            const double eri_ijkl = Lambda * (

                    S_ij_01 * S_kl_00 * F1_t[0] +

                    S_ij_00 * S_kl_00 * ((a_k + a_l) / (a_i + a_j + a_k + a_l)) * rPQ[j_cart] * F1_t[1]

                    );

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * ss_mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < sp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockSPSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   sp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < sp_prim_pair_count) && (ij < sp_prim_pair_count) && (std::fabs(sp_mat_Q[ij] * sp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto k = sp_first_inds[kl];
            const auto l = sp_second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto j_cart = j % 3;
            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 2, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};
            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;
            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_ij_01 * S_kl_01 * F2_t[0] +

                    S_ij_01 * S_kl_00 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[l_cart] * F2_t[1] +

                    S_ij_00 * S_kl_01 * ((a_k + a_l) / (a_i + a_j + a_k + a_l)) * rPQ[j_cart] * F2_t[1] +

                    S_ij_00 * S_kl_00 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[l_cart] * F2_t[2] +

                        (j_cart == l_cart ? 1.0 : 0.0) * 0.5 * F2_t[1]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * sp_mat_D[kl] * 2.0;
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < sp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockSPPP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   sp_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < pp_prim_pair_count) && (ij < sp_prim_pair_count) && (std::fabs(sp_mat_Q[ij] * pp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto k = pp_first_inds[kl];
            const auto l = pp_second_inds[kl];

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

            const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];
            const auto c_k = p_prim_info[k / 3 + p_prim_count * 1];
            const auto x_k = p_prim_info[k / 3 + p_prim_count * 2];
            const auto y_k = p_prim_info[k / 3 + p_prim_count * 3];
            const auto z_k = p_prim_info[k / 3 + p_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto j_cart = j % 3;
            const auto k_cart = k % 3;
            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 3, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};
            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

            const auto S_kl_10 = (a_l / (a_k + a_l)) * rkl[k_cart] * S_kl_00;
            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const auto S_kl_11 = (

                    (k_cart == l_cart ? 1.0 : 0.0) * (0.5 / (a_k + a_l)) +

                    (a_k * a_l / ((a_k + a_l) * (a_k + a_l))) * (-rkl[k_cart]) * rkl[l_cart]

                    ) * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_ij_01 * S_kl_11 * F3_t[0] +

                    S_ij_00 * S_kl_11 * ((a_k + a_l) / (a_i + a_j + a_k + a_l)) * rPQ[j_cart] * F3_t[1] +

                    S_ij_01 * S_kl_01 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[k_cart] * F3_t[1] +

                    S_ij_00 * S_kl_01 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[k_cart] * F3_t[2] +

                        (j_cart == k_cart ? 1.0 : 0.0) * 0.5 * F3_t[1]

                    ) +

                    S_ij_01 * S_kl_10 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[l_cart] * F3_t[1] +

                    S_ij_00 * S_kl_10 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[l_cart] * F3_t[2] +

                        (j_cart == l_cart ? 1.0 : 0.0) * 0.5 * F3_t[1]

                    ) +

                    S_ij_01 * S_kl_00 * (a_i + a_j) / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) / (a_i + a_j + a_k + a_l) * rPQ[k_cart] * rPQ[l_cart] * F3_t[2] -

                        (k_cart == l_cart ? 1.0 : 0.0) * 0.5 / (a_k + a_l) * F3_t[1]

                    ) +

                    S_ij_00 * S_kl_00 * (-(a_i + a_j)) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[k_cart] * rPQ[l_cart] * F3_t[3] -

                        ((j_cart == k_cart ? 1.0 : 0.0) * (-rPQ[l_cart]) + 
                         (j_cart == l_cart ? 1.0 : 0.0) * (-rPQ[k_cart]) +
                         (k_cart == l_cart ? 1.0 : 0.0) * (-rPQ[j_cart])) * 0.5 * F3_t[2]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * pp_mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < sp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockPPSS(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   ss_mat_D,
                       const double*   ss_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < ss_prim_pair_count) && (ij < pp_prim_pair_count) && (std::fabs(pp_mat_Q[ij] * ss_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            const auto k = ss_first_inds[kl];
            const auto l = ss_second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = s_prim_info[l + s_prim_count * 0];
            const auto c_l = s_prim_info[l + s_prim_count * 1];
            const auto x_l = s_prim_info[l + s_prim_count * 2];
            const auto y_l = s_prim_info[l + s_prim_count * 3];
            const auto z_l = s_prim_info[l + s_prim_count * 4];

            const auto i_cart = i % 3;
            const auto j_cart = j % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 2, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};

            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;

            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

            const auto S_ij_11 = (

                    (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) +

                    (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                    ) * S_ij_00;

            const double eri_ijkl = Lambda * (

                    S_kl_00 * S_ij_11 * F2_t[0] +

                    S_kl_00 * S_ij_10 * ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * rPQ[j_cart] * F2_t[1] +

                    S_kl_00 * S_ij_01 * ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * rPQ[i_cart] * F2_t[1] +

                    S_kl_00 * S_ij_00 * ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * (

                        ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * rPQ[i_cart] * rPQ[j_cart] * F2_t[2] -

                        (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) * F2_t[1]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * ss_mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < pp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockPPSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   sp_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (sp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < sp_prim_pair_count) && (ij < pp_prim_pair_count) && (std::fabs(pp_mat_Q[ij] * sp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            const auto k = sp_first_inds[kl];
            const auto l = sp_second_inds[kl];

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

            const auto a_k = s_prim_info[k + s_prim_count * 0];
            const auto c_k = s_prim_info[k + s_prim_count * 1];
            const auto x_k = s_prim_info[k + s_prim_count * 2];
            const auto y_k = s_prim_info[k + s_prim_count * 3];
            const auto z_k = s_prim_info[k + s_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto i_cart = i % 3;
            const auto j_cart = j % 3;
            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 3, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};
            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;
            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

            const auto S_ij_11 = (

                    (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) +

                    (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                    ) * S_ij_00;

            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_kl_01 * S_ij_11 * F3_t[0] +

                    S_kl_00 * S_ij_11 * (-(a_i + a_j) / (a_k + a_l + a_i + a_j)) * rPQ[l_cart] * F3_t[1] +

                    S_kl_01 * S_ij_10 * ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * rPQ[j_cart] * F3_t[1] +

                    S_kl_00 * S_ij_10 / (a_k + a_l + a_i + a_j) * (

                        (a_k + a_l) * (a_i + a_j) / (a_k + a_l + a_i + a_j) * (-rPQ[j_cart]) * rPQ[l_cart] * F3_t[2] +

                        (j_cart == l_cart ? 1.0 : 0.0) * 0.5 * F3_t[1]

                    ) +

                    S_kl_01 * S_ij_01 * ((a_k + a_l) / (a_k + a_l + a_i + a_j)) * rPQ[i_cart] * F3_t[1] +

                    S_kl_00 * S_ij_01 / (a_k + a_l + a_i + a_j) * (

                        (a_k + a_l) * (a_i + a_j) / (a_k + a_l + a_i + a_j) * (-rPQ[i_cart]) * rPQ[l_cart] * F3_t[2] +

                        (i_cart == l_cart ? 1.0 : 0.0) * 0.5 * F3_t[1]

                    ) +

                    S_kl_01 * S_ij_00 * (a_k + a_l) / (a_k + a_l + a_i + a_j) * (

                        (a_k + a_l) / (a_k + a_l + a_i + a_j) * rPQ[i_cart] * rPQ[j_cart] * F3_t[2] -

                        (i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F3_t[1]

                    ) +

                    S_kl_00 * S_ij_00 * (-(a_k + a_l)) / ((a_k + a_l + a_i + a_j) * (a_k + a_l + a_i + a_j)) * (

                        (a_k + a_l) * (a_i + a_j) / (a_k + a_l + a_i + a_j) * rPQ[i_cart] * rPQ[j_cart] * rPQ[l_cart] * F3_t[3] +

                        ((l_cart == i_cart ? 1.0 : 0.0) * (-rPQ[j_cart]) + 
                         (l_cart == j_cart ? 1.0 : 0.0) * (-rPQ[i_cart]) +
                         (i_cart == j_cart ? 1.0 : 0.0) * (-rPQ[l_cart])) * 0.5 * F3_t[2]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * sp_mat_D[kl] * 2.0;
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < pp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
computeCoulombFockPPPP(double*         mat_J,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   pp_mat_Q,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double J_ij = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

    __syncthreads();

    for (uint32_t m = 0; m < (pp_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < pp_prim_pair_count) && (ij < pp_prim_pair_count) && (std::fabs(pp_mat_Q[ij] * pp_mat_Q[kl] * max_D) > 1.0e-12))
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            const auto k = pp_first_inds[kl];
            const auto l = pp_second_inds[kl];

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

            const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];
            const auto c_k = p_prim_info[k / 3 + p_prim_count * 1];
            const auto x_k = p_prim_info[k / 3 + p_prim_count * 2];
            const auto y_k = p_prim_info[k / 3 + p_prim_count * 3];
            const auto z_k = p_prim_info[k / 3 + p_prim_count * 4];

            const auto a_l = p_prim_info[l / 3 + p_prim_count * 0];
            const auto c_l = p_prim_info[l / 3 + p_prim_count * 1];
            const auto x_l = p_prim_info[l / 3 + p_prim_count * 2];
            const auto y_l = p_prim_info[l / 3 + p_prim_count * 3];
            const auto z_l = p_prim_info[l / 3 + p_prim_count * 4];

            const auto i_cart = i % 3;
            const auto j_cart = j % 3;
            const auto k_cart = k % 3;
            const auto l_cart = l % 3;

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
            const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rPQ[3] = {(a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j),
                                   (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j),
                                   (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j)};

            const auto r2_PQ = rPQ[0] * rPQ[0] + rPQ[1] * rPQ[1] + rPQ[2] * rPQ[2];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

            double F4_t[5];

            gpu::computeBoysFunction(F4_t, (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ, 4, boys_func_table, boys_func_ft);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl_00 = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};
            const double rkl[3] = {y_l - y_k, z_l - z_k, x_l - x_k};

            const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;
            const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

            const auto S_ij_11 = (

                    (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) +

                    (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                    ) * S_ij_00;

            const auto S_kl_10 = (a_l / (a_k + a_l)) * rkl[k_cart] * S_kl_00;
            const auto S_kl_01 = -(a_k / (a_k + a_l)) * rkl[l_cart] * S_kl_00;

            const auto S_kl_11 = (

                    (k_cart == l_cart ? 1.0 : 0.0) * (0.5 / (a_k + a_l)) +

                    (a_k * a_l / ((a_k + a_l) * (a_k + a_l))) * (-rkl[k_cart]) * rkl[l_cart]

                    ) * S_kl_00;

            const double eri_ijkl = Lambda * (

                    S_ij_11 * S_kl_11 * F4_t[0]

                    + S_ij_11 * S_kl_10 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[l_cart] * F4_t[1] +

                    + S_ij_11 * S_kl_01 * (-(a_i + a_j) / (a_i + a_j + a_k + a_l)) * rPQ[k_cart] * F4_t[1] +

                    + S_ij_11 * S_kl_00 * (a_i + a_j) / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) / (a_i + a_j + a_k + a_l) * rPQ[k_cart] * rPQ[l_cart] * F4_t[2]

                        - (k_cart == l_cart ? 1.0 : 0.0) * 0.5 / (a_k + a_l) * F4_t[1]

                    )

                    + S_ij_10 * S_kl_11 * ((a_k + a_l) / (a_i + a_j + a_k + a_l)) * rPQ[j_cart] * F4_t[1] +

                    + S_ij_10 * S_kl_10 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[l_cart] * F4_t[2]

                        + (j_cart == l_cart ? 1.0 : 0.0) * 0.5 * F4_t[1]

                    )

                    + S_ij_10 * S_kl_01 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[k_cart] * F4_t[2]

                        + (j_cart == k_cart ? 1.0 : 0.0) * 0.5 * F4_t[1]

                    )

                    + S_ij_10 * S_kl_00 * (-(a_i + a_j)) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[j_cart]) * rPQ[k_cart] * rPQ[l_cart] * F4_t[3]

                        - ((j_cart == k_cart ? 1.0 : 0.0) * (-rPQ[l_cart]) + 
                           (j_cart == l_cart ? 1.0 : 0.0) * (-rPQ[k_cart]) +
                           (k_cart == l_cart ? 1.0 : 0.0) * (-rPQ[j_cart])) * 0.5 * F4_t[2]

                    )

                    + S_ij_01 * S_kl_11 * ((a_k + a_l) / (a_i + a_j + a_k + a_l)) * rPQ[i_cart] * F4_t[1] +

                    + S_ij_01 * S_kl_10 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[i_cart]) * rPQ[l_cart] * F4_t[2]

                        + (i_cart == l_cart ? 1.0 : 0.0) * 0.5 * F4_t[1]

                    )

                    + S_ij_01 * S_kl_01 / (a_i + a_j + a_k + a_l) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[i_cart]) * rPQ[k_cart] * F4_t[2]

                        + (i_cart == k_cart ? 1.0 : 0.0) * 0.5 * F4_t[1]

                    )

                    + S_ij_01 * S_kl_00 * (-(a_i + a_j)) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (-rPQ[i_cart]) * rPQ[k_cart] * rPQ[l_cart] * F4_t[3]

                        - ((i_cart == k_cart ? 1.0 : 0.0) * (-rPQ[l_cart]) + 
                           (i_cart == l_cart ? 1.0 : 0.0) * (-rPQ[k_cart]) +
                           (k_cart == l_cart ? 1.0 : 0.0) * (-rPQ[i_cart])) * 0.5 * F4_t[2]

                    )

                    + S_ij_00 * S_kl_11 * (a_k + a_l) / (a_i + a_j + a_k + a_l) * (

                        (a_k + a_l) / (a_i + a_j + a_k + a_l) * rPQ[i_cart] * rPQ[j_cart] * F4_t[2]

                        - (i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F4_t[1]

                    )

                    + S_ij_00 * S_kl_10 * (-(a_k + a_l)) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * rPQ[i_cart] * rPQ[j_cart] * rPQ[l_cart] * F4_t[3]

                        + ((i_cart == j_cart ? 1.0 : 0.0) * (-rPQ[l_cart]) + 
                           (i_cart == l_cart ? 1.0 : 0.0) * (-rPQ[j_cart]) +
                           (j_cart == l_cart ? 1.0 : 0.0) * (-rPQ[i_cart])) * 0.5 * F4_t[2]

                    )

                    + S_ij_00 * S_kl_01 * (-(a_k + a_l)) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * rPQ[i_cart] * rPQ[j_cart] * rPQ[k_cart] * F4_t[3]

                        + ((i_cart == j_cart ? 1.0 : 0.0) * (-rPQ[k_cart]) + 
                           (i_cart == k_cart ? 1.0 : 0.0) * (-rPQ[j_cart]) +
                           (j_cart == k_cart ? 1.0 : 0.0) * (-rPQ[i_cart])) * 0.5 * F4_t[2]

                    )

                    + S_ij_00 * S_kl_00 / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * (

                        (a_i + a_j) * (a_i + a_j) * (a_k + a_l) * (a_k + a_l) / ((a_i + a_j + a_k + a_l) * (a_i + a_j + a_k + a_l)) * rPQ[i_cart] * rPQ[j_cart] * rPQ[k_cart] * rPQ[l_cart] * F4_t[4]

                        - (a_i + a_j) * (a_k + a_l) / ((a_i + a_j + a_k + a_l) * 2.0) * (

                              (i_cart == j_cart ? 1.0 : 0.0) * rPQ[k_cart] * rPQ[l_cart]
                            + (i_cart == k_cart ? 1.0 : 0.0) * rPQ[j_cart] * rPQ[l_cart]
                            + (i_cart == l_cart ? 1.0 : 0.0) * rPQ[j_cart] * rPQ[k_cart]
                            + (j_cart == k_cart ? 1.0 : 0.0) * rPQ[i_cart] * rPQ[l_cart]
                            + (j_cart == l_cart ? 1.0 : 0.0) * rPQ[i_cart] * rPQ[k_cart]
                            + (k_cart == l_cart ? 1.0 : 0.0) * rPQ[i_cart] * rPQ[j_cart]

                        ) * F4_t[3]

                        + ((i_cart == j_cart ? 1.0 : 0.0) * (k_cart == l_cart ? 1.0 : 0.0) + 
                           (i_cart == k_cart ? 1.0 : 0.0) * (j_cart == l_cart ? 1.0 : 0.0) +
                           (i_cart == l_cart ? 1.0 : 0.0) * (j_cart == k_cart ? 1.0 : 0.0)) * 0.25 * F4_t[2]

                    ));

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * pp_mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;

            // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
            if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
        }

        __syncthreads();

        // early exit for thread block
        if (skip_thread_block == 1) break;

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();
    }

    if ((threadIdx.y == 0) && (ij < pp_prim_pair_count))
    {
        mat_J[ij] = J_ij;
    }
}

__global__ void
cudaFockK(double*         mat_K,
          const uint32_t* pair_inds_i,
          const uint32_t* pair_inds_k,
          const double*   s_prim_info,
          const uint32_t* s_prim_aoinds,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D_full,
          const double*   mat_Q_for_K,
          const uint32_t* first_inds_for_K,
          const uint32_t* second_inds_for_K,
          const uint32_t  pair_inds_count,
          const uint32_t  naos,
          const double*   boys_func_table,
          const double*   boys_func_ft)
{
    // each thread block scans over [i?|k?] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ik = blockIdx.x;

    double K_ik = 0.0;

    for (uint32_t m = 0; m < (s_prim_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t j = m * TILE_DIM + threadIdx.y;

        // to avoid memory hazard
        __syncthreads();

        if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

        __syncthreads();

        for (uint32_t n = 0; n < (s_prim_count + TILE_DIM - 1) / TILE_DIM; n++)
        {
            const uint32_t l = n * TILE_DIM + threadIdx.x;

            if ((ik < pair_inds_count) && (j < s_prim_count) && (l < s_prim_count))
            {
                const auto i = pair_inds_i[ik];
                const auto k = pair_inds_k[ik];

                const auto Q_ij = mat_Q_for_K[i * s_prim_count + j];
                const auto Q_kl = mat_Q_for_K[k * s_prim_count + l];

                if (std::fabs(Q_ij * Q_kl * max_D) > 1.0e-12)
                {
                    auto j_prim = second_inds_for_K[i * s_prim_count + j];
                    auto l_prim = second_inds_for_K[k * s_prim_count + l];

                    auto j_cgto = s_prim_aoinds[j_prim];
                    auto l_cgto = s_prim_aoinds[l_prim];

                    const auto a_i = s_prim_info[i + s_prim_count * 0];
                    const auto c_i = s_prim_info[i + s_prim_count * 1];
                    const auto x_i = s_prim_info[i + s_prim_count * 2];
                    const auto y_i = s_prim_info[i + s_prim_count * 3];
                    const auto z_i = s_prim_info[i + s_prim_count * 4];

                    const auto a_j = s_prim_info[j_prim + s_prim_count * 0];
                    const auto c_j = s_prim_info[j_prim + s_prim_count * 1];
                    const auto x_j = s_prim_info[j_prim + s_prim_count * 2];
                    const auto y_j = s_prim_info[j_prim + s_prim_count * 3];
                    const auto z_j = s_prim_info[j_prim + s_prim_count * 4];

                    const auto a_k = s_prim_info[k + s_prim_count * 0];
                    const auto c_k = s_prim_info[k + s_prim_count * 1];
                    const auto x_k = s_prim_info[k + s_prim_count * 2];
                    const auto y_k = s_prim_info[k + s_prim_count * 3];
                    const auto z_k = s_prim_info[k + s_prim_count * 4];

                    const auto a_l = s_prim_info[l_prim + s_prim_count * 0];
                    const auto c_l = s_prim_info[l_prim + s_prim_count * 1];
                    const auto x_l = s_prim_info[l_prim + s_prim_count * 2];
                    const auto y_l = s_prim_info[l_prim + s_prim_count * 3];
                    const auto z_l = s_prim_info[l_prim + s_prim_count * 4];

                    const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);
                    const auto r2_kl = (x_l - x_k) * (x_l - x_k) + (y_l - y_k) * (y_l - y_k) + (z_l - z_k) * (z_l - z_k);

                    const auto x_PQ = (a_k * x_k + a_l * x_l) / (a_k + a_l) - (a_i * x_i + a_j * x_j) / (a_i + a_j);
                    const auto y_PQ = (a_k * y_k + a_l * y_l) / (a_k + a_l) - (a_i * y_i + a_j * y_j) / (a_i + a_j);
                    const auto z_PQ = (a_k * z_k + a_l * z_l) / (a_k + a_l) - (a_i * z_i + a_j * z_j) / (a_i + a_j);

                    const auto r2_PQ = x_PQ * x_PQ + y_PQ * y_PQ + z_PQ * z_PQ;

                    // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

                    const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (MATH_CONST_PI * (a_i + a_j + a_k + a_l)));

                    const auto t = (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ;

                    double F0_t[1];

                    gpu::computeBoysFunction(F0_t, t, 0, boys_func_table, boys_func_ft);

                    const auto S_ij = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
                    const auto S_kl = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

                    ERIs[threadIdx.y][threadIdx.x] = Lambda * S_ij * S_kl * F0_t[0] * mat_D_full[j_cgto * naos + l_cgto];
                }
                else
                {
                    ERIs[threadIdx.y][threadIdx.x] = 0.0;

                    // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
                    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = m + 1;
                }
            }
            else
            {
                ERIs[threadIdx.y][threadIdx.x] = 0.0;

                // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
                if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = m + 1;
            }

            __syncthreads();

            // early exit for thread block
            if (skip_thread_block == m + 1) break;

            if ((threadIdx.y == 0) && (threadIdx.x == 0))
            {
                for (uint32_t y = 0; y < TILE_DIM; y++)
                {
                    for (uint32_t x = 0; x < TILE_DIM; x++)
                    {
                        K_ik += ERIs[y][x];
                    }
                }
            }

            __syncthreads();
        }
    }

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count))
    {
        mat_K[ik] = K_ik;
    }
}

}  // namespace gpu
