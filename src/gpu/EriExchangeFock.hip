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
#include "EriExchangeFock.hpp"
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
cudaFockK(double*         mat_K,
          const uint32_t* pari_inds_i,
          const uint32_t* pari_inds_k,
          const double*   s_prim_info,
          const uint32_t* s_prim_aoinds,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D_full,
          const double*   mat_Q_for_K,
          const uint32_t* first_inds_for_K,
          const uint32_t* second_inds_for_K,
          const uint32_t  ss_prim_pair_count,
          const uint32_t  naos)
{
    // each thread block scans over [i?|k?] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t skip_thread_block;

    const uint32_t ik = blockIdx.x;

    double K_ik = 0.0;

    for (uint32_t m = 0; m < (s_prim_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t j = m * TILE_DIM + threadIdx.y;

        skip_thread_block = 0;

        for (uint32_t n = 0; n < (s_prim_count + TILE_DIM - 1) / TILE_DIM; n++)
        {
            const uint32_t l = n * TILE_DIM + threadIdx.x;

            if ((ik < ss_prim_pair_count) && (j < s_prim_count) && (l < s_prim_count))
            {
                const auto i = pari_inds_i[ik];
                const auto k = pari_inds_k[ik];

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

                    // TODO: higher order Boys function

                    const double F0_t = (t == 0.0 ? 1.0 : std::sqrt(MATH_CONST_PI / (4.0 * t)) * std::erf(std::sqrt(t)));

                    const auto S_ij = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
                    const auto S_kl = c_k * c_l * std::pow(MATH_CONST_PI / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

                    ERIs[threadIdx.y][threadIdx.x] = Lambda * S_ij * S_kl * F0_t * mat_D_full[j_cgto * naos + l_cgto];
                }
                else
                {
                    ERIs[threadIdx.y][threadIdx.x] = 0.0;

                    // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
                    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 1;
                }
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

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < ss_prim_pair_count))
    {
        mat_K[ik] = K_ik;
    }
}

}  // namespace gpu
