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

#include "EriSSSS.hpp"
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

#define hipSafe(e)                                                                                                       \
    {                                                                                                                     \
        hipError_t err = (e);                                                                                            \
        if (err != hipSuccess)                                                                                           \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                      \
        }                                                                                                                 \
    }

namespace gpu {  // gpu namespace

__global__ void
cudaFockJ(double*         mat_J,
          const double*   s_prim_info,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D,
          const double*   mat_Q,
          const uint32_t* first_inds,
          const uint32_t* second_inds,
          const uint32_t  ss_prim_pair_count)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    const double m_pi = 3.14159265358979323846;

    double J_ij = 0.0;

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        // TODO: add D_kl into ket side

        if ((kl < ss_prim_pair_count) && (ij < ss_prim_pair_count) && (mat_Q[ij] * mat_Q[kl] * max_D > 1.0e-15))
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

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (m_pi * (a_i + a_j + a_k + a_l)));

            const auto t = (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ;

            // TODO: higher order Boys function

            const double F0_t = (t == 0.0 ? 1.0 : std::sqrt(m_pi / (4.0 * t)) * std::erf(std::sqrt(t)));

            const auto S_ij = c_i * c_j * std::pow(m_pi / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
            const auto S_kl = c_k * c_l * std::pow(m_pi / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

            // NOTE: doubling for off-diagonal elements of D due to k<=>l symmetry
            ERIs[threadIdx.y][threadIdx.x] = Lambda * S_ij * S_kl * F0_t * mat_D[kl] * (k == l ? 1.0 : 2.0);
        }
        else
        {
            ERIs[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

        // early exit for thread block (ERIs[0][0] has the largest upper bound)
        if (ERIs[0][0] == 0.0) break;

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

    const uint32_t ik = blockIdx.x;

    const double m_pi = 3.14159265358979323846;

    double K_ik = 0.0;

    for (uint32_t m = 0; m < (s_prim_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t j = m * TILE_DIM + threadIdx.y;

        for (uint32_t n = 0; n < (s_prim_count + TILE_DIM - 1) / TILE_DIM; n++)
        {
            const uint32_t l = n * TILE_DIM + threadIdx.x;

            if ((ik < ss_prim_pair_count) && (j < s_prim_count) && (l < s_prim_count))
            {
                const auto i = pari_inds_i[ik];
                const auto k = pari_inds_k[ik];

                const auto Q_ij = mat_Q_for_K[i * s_prim_count + j];
                const auto Q_kl = mat_Q_for_K[k * s_prim_count + l];

                if (Q_ij * Q_kl * max_D > 1.0e-15)
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

                    const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_k + a_l) / (m_pi * (a_i + a_j + a_k + a_l)));

                    const auto t = (a_i + a_j) * (a_k + a_l) / (a_i + a_j + a_k + a_l) * r2_PQ;

                    // TODO: higher order Boys function

                    const double F0_t = (t == 0.0 ? 1.0 : std::sqrt(m_pi / (4.0 * t)) * std::erf(std::sqrt(t)));

                    const auto S_ij = c_i * c_j * std::pow(m_pi / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);
                    const auto S_kl = c_k * c_l * std::pow(m_pi / (a_k + a_l), 1.5) * std::exp(-a_k * a_l / (a_k + a_l) * r2_kl);

                    ERIs[threadIdx.y][threadIdx.x] = Lambda * S_ij * S_kl * F0_t * mat_D_full[j_cgto * naos + l_cgto];
                }
                else
                {
                    ERIs[threadIdx.y][threadIdx.x] = 0.0;
                }
            }
            else
            {
                ERIs[threadIdx.y][threadIdx.x] = 0.0;
            }

            __syncthreads();

            // early exit for thread block (ERIs[0][0] has the largest upper bound)
            if (ERIs[0][0] == 0.0) break;

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

auto
computeEriSSSS(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> CDenseMatrix
{
    const auto m_pi = mathconst::getPiValue();

    auto rank = mpi::rank(MPI_COMM_WORLD);

    if (rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (rank % 8 == 7) hipSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::computeEriSSSS: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    // S gto block

    int64_t s_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        if (gto_block.getAngularMomentum() == 0)
        {
            const auto ncgtos = gto_block.getNumberOfBasisFunctions();
            const auto npgtos = gto_block.getNumberOfPrimitives();

            s_prim_count += npgtos * ncgtos;
        }
    }

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    int64_t s_prim_offset = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0)
        {
            const auto ncgtos = gto_block.getNumberOfBasisFunctions();
            const auto npgtos = gto_block.getNumberOfPrimitives();

            const auto gto_exps   = gto_block.getExponents();
            const auto gto_norms  = gto_block.getNormalizationFactors();
            const auto gto_coords = gto_block.getCoordinates();

            const auto gto_ao_inds = gto_block.getAtomicOrbitalsIndexes();

            for (int64_t j = 0; j < npgtos; j++)
            {
                for (int64_t i = 0; i < ncgtos; i++)
                {
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = gto_exps[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 1] = gto_norms[i + j * ncgtos];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 2] = gto_coords[i][0];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 3] = gto_coords[i][1];
                    s_prim_info[s_prim_offset + i + j * ncgtos + s_prim_count * 4] = gto_coords[i][2];

                    s_prim_aoinds[s_prim_offset + i + j * ncgtos + s_prim_count * 0] = static_cast<uint32_t>(gto_ao_inds[i]);
                }
            }

            s_prim_offset += npgtos * ncgtos;
        }
    }

    double*   d_s_prim_info;
    uint32_t* d_s_prim_aoinds;

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_s_prim_aoinds, s_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // For now, only process S gto block

    // TODO: higher angular momentum

    auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;

    CDenseMatrix mat_Q_orig(s_prim_count, s_prim_count);

    std::vector<std::tuple<double, int64_t, int64_t, double>> mat_Q_D_sorted;

    auto dens_ptr = densityMatrix.alphaDensity(0);

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = i; j < s_prim_count; j++, ij++)
        {
            const auto a_j = s_prim_info[j + s_prim_count * 0];
            const auto c_j = s_prim_info[j + s_prim_count * 1];
            const auto x_j = s_prim_info[j + s_prim_count * 2];
            const auto y_j = s_prim_info[j + s_prim_count * 3];
            const auto z_j = s_prim_info[j + s_prim_count * 4];

            const auto j_cgto = s_prim_aoinds[j];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij = c_i * c_j * std::pow(m_pi / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            // Note: t = 0.0 and F0_t = 1.0
            const auto eri_ijij = S_ij * S_ij * std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / m_pi);

            const auto sqrt_ijij = std::sqrt(eri_ijij);

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

            mat_Q_orig.row(i)[j] = sqrt_ijij;

            if (i != j) mat_Q_orig.row(j)[i] = sqrt_ijij;

            mat_Q_D_sorted.push_back(std::make_tuple(sqrt_ijij, i, j, D_ij));
        }
    }

    std::sort(mat_Q_D_sorted.begin(), mat_Q_D_sorted.end());
    std::reverse(mat_Q_D_sorted.begin(), mat_Q_D_sorted.end());

    // K preparation

    std::vector<double>   mat_Q_for_K(s_prim_count * s_prim_count);
    std::vector<uint32_t> first_inds_for_K(s_prim_count * s_prim_count);
    std::vector<uint32_t> second_inds_for_K(s_prim_count * s_prim_count);

    // TODO: use uplo_index

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = mat_Q_orig.row(i);

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            mat_Q_for_K[i * s_prim_count + j]       = q_val;
            first_inds_for_K[i * s_prim_count + j]  = static_cast<uint32_t>(i_idx);
            second_inds_for_K[i * s_prim_count + j] = static_cast<uint32_t>(j_idx);
        }
    }

    std::vector<uint32_t> pair_inds_i(ss_prim_pair_count);
    std::vector<uint32_t> pair_inds_k(ss_prim_pair_count);

    for (int64_t i = 0, ik = 0; i < s_prim_count; i++)
    {
        for (int64_t k = i; k < s_prim_count; k++, ik++)
        {
            pair_inds_i[ik] = i;
            pair_inds_k[ik] = k;
        }
    }

    // sorted Q, D, and indices on host

    std::vector<double> mat_Q(ss_prim_pair_count);
    std::vector<double> mat_D(ss_prim_pair_count);
    std::vector<double> mat_J(ss_prim_pair_count);

    std::vector<uint32_t> first_inds(ss_prim_pair_count);
    std::vector<uint32_t> second_inds(ss_prim_pair_count);

    for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto& q_ij = mat_Q_D_sorted[ij];

        auto sqrt_ijij = std::get<0>(q_ij);
        auto i         = std::get<1>(q_ij);
        auto j         = std::get<2>(q_ij);
        auto D_ij      = std::get<3>(q_ij);

        mat_Q[ij] = sqrt_ijij;
        mat_D[ij] = D_ij;

        first_inds[ij]  = static_cast<uint32_t>(i);
        second_inds[ij] = static_cast<uint32_t>(j);
    }

    auto max_D_ptr = std::max_element(mat_D.begin(), mat_D.end());

    // sorted Q, D, and indices on device

    double *  d_mat_J, *d_mat_D, *d_mat_Q;
    uint32_t *d_first_inds, *d_second_inds;

    hipSafe(hipMalloc(&d_mat_Q, ss_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_D, ss_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_J, ss_prim_pair_count * sizeof(double)));

    hipSafe(hipMalloc(&d_first_inds, ss_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_second_inds, ss_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_mat_D, mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q, mat_Q.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_first_inds, first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_second_inds, second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    double*   d_mat_K;
    uint32_t *d_pair_inds_i, *d_pair_inds_k;

    hipSafe(hipMalloc(&d_mat_K, ss_prim_pair_count * sizeof(double)));

    hipSafe(hipMalloc(&d_pair_inds_i, ss_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k, ss_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_pair_inds_i, pair_inds_i.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k, pair_inds_k.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    double *  d_mat_Q_for_K, *d_mat_D_full;
    uint32_t *d_first_inds_for_K, *d_second_inds_for_K;

    hipSafe(hipMalloc(&d_mat_D_full, naos * naos * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Q_for_K, s_prim_count * s_prim_count * sizeof(double)));

    hipSafe(hipMalloc(&d_first_inds_for_K, s_prim_count * s_prim_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_second_inds_for_K, s_prim_count * s_prim_count * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_mat_D_full, dens_ptr, naos * naos * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q_for_K, mat_Q_for_K.data(), s_prim_count * s_prim_count * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_first_inds_for_K, first_inds_for_K.data(), s_prim_count * s_prim_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_second_inds_for_K, second_inds_for_K.data(), s_prim_count * s_prim_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("Preparation");

    timer.start("J computation");

    // compute J

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((ss_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

    hipLaunchKernelGGL(cudaFockJ, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                 d_s_prim_info,
                                                 static_cast<uint32_t>(s_prim_count),
                                                 *max_D_ptr,
                                                 d_mat_D,
                                                 d_mat_Q,
                                                 d_first_inds,
                                                 d_second_inds,
                                                 static_cast<uint32_t>(ss_prim_pair_count));

    hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

    timer.stop("J computation");

    timer.start("K computation");

    threads_per_block = dim3(TILE_DIM, TILE_DIM);

    num_blocks = dim3(ss_prim_pair_count, 1);

    hipLaunchKernelGGL(cudaFockK, num_blocks, threads_per_block, 0, 0, d_mat_K,
                                                 d_pair_inds_i,
                                                 d_pair_inds_k,
                                                 d_s_prim_info,
                                                 d_s_prim_aoinds,
                                                 static_cast<uint32_t>(s_prim_count),
                                                 *max_D_ptr,
                                                 d_mat_D_full,
                                                 d_mat_Q_for_K,
                                                 d_first_inds_for_K,
                                                 d_second_inds_for_K,
                                                 static_cast<uint32_t>(ss_prim_pair_count),
                                                 static_cast<uint32_t>(naos));

    std::vector<double> mat_K(ss_prim_pair_count);

    hipSafe(hipMemcpy(mat_K.data(), d_mat_K, ss_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

    timer.stop("K computation");

    timer.start("J/K to Fock");

    // put J back to Fock

    CDenseMatrix mat_Fock(naos, naos);

    mat_Fock.zero();

    for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto i = first_inds[ij];
        const auto j = second_inds[ij];

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = s_prim_aoinds[j];

        mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

        if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
    }

    for (int64_t ik = 0; ik < ss_prim_pair_count; ik++)
    {
        const auto i = pair_inds_i[ik];
        const auto k = pair_inds_k[ik];

        const auto i_cgto = s_prim_aoinds[i];
        const auto k_cgto = s_prim_aoinds[k];

        mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

        if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
    }

    timer.stop("J/K to Fock");

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_s_prim_aoinds));

    hipSafe(hipFree(d_mat_Q));
    hipSafe(hipFree(d_mat_D));
    hipSafe(hipFree(d_mat_J));

    hipSafe(hipFree(d_first_inds));
    hipSafe(hipFree(d_second_inds));

    hipSafe(hipFree(d_mat_K));

    hipSafe(hipFree(d_pair_inds_i));
    hipSafe(hipFree(d_pair_inds_k));

    hipSafe(hipFree(d_mat_D_full));
    hipSafe(hipFree(d_mat_Q_for_K));

    hipSafe(hipFree(d_first_inds_for_K));
    hipSafe(hipFree(d_second_inds_for_K));

    timer.stop("Total timing");

    std::cout << "\nTiming of ERIs on GPU\n";
    std::cout << "-----------------------\n";
    std::cout << timer.getSummary() << std::endl;

    return mat_Fock;
}

}  // namespace gpu
