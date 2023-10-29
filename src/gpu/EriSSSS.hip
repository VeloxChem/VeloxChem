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

#define cudaSafe(e)                                                                                                       \
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
cuda_mat_J_PI(double*         mat_J,
              const double*   gto_info,
              const uint32_t  n_primitives,
              const double*   mat_D,
              const double*   mat_Q,
              const uint32_t* first_inds,
              const uint32_t* second_inds,
              const uint32_t  n_prim_pairs)
{
    // each thread block scans over [ij|??] and sum up to a primitive J matrix element

    __shared__ double ERIs[TILE_DIM][TILE_DIM + 1];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    const double m_pi = 3.14159265358979323846;

    double J_ij = 0.0;

    for (uint32_t m = 0; m < (n_prim_pairs + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl < n_prim_pairs) && (ij < n_prim_pairs) && (mat_Q[ij] * mat_Q[kl] > 1.0e-15))
        {
            const auto i = first_inds[ij];
            const auto j = second_inds[ij];

            const auto k = first_inds[kl];
            const auto l = second_inds[kl];

            const auto a_i = gto_info[i + n_primitives * 0];
            const auto c_i = gto_info[i + n_primitives * 1];
            const auto x_i = gto_info[i + n_primitives * 2];
            const auto y_i = gto_info[i + n_primitives * 3];
            const auto z_i = gto_info[i + n_primitives * 4];

            const auto a_j = gto_info[j + n_primitives * 0];
            const auto c_j = gto_info[j + n_primitives * 1];
            const auto x_j = gto_info[j + n_primitives * 2];
            const auto y_j = gto_info[j + n_primitives * 3];
            const auto z_j = gto_info[j + n_primitives * 4];

            const auto a_k = gto_info[k + n_primitives * 0];
            const auto c_k = gto_info[k + n_primitives * 1];
            const auto x_k = gto_info[k + n_primitives * 2];
            const auto y_k = gto_info[k + n_primitives * 3];
            const auto z_k = gto_info[k + n_primitives * 4];

            const auto a_l = gto_info[l + n_primitives * 0];
            const auto c_l = gto_info[l + n_primitives * 1];
            const auto x_l = gto_info[l + n_primitives * 2];
            const auto y_l = gto_info[l + n_primitives * 3];
            const auto z_l = gto_info[l + n_primitives * 4];

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

        if (threadIdx.y == 0)
        {
            for (uint32_t n = 0; n < TILE_DIM; n++)
            {
                J_ij += ERIs[n][threadIdx.x];
            }
        }

        __syncthreads();

        // Note: early exit for thread block (ERIs[0][0] has the largest upper bound)
        if (ERIs[0][0] == 0.0) break;
    }

    if ((threadIdx.y == 0) && (ij < n_prim_pairs))
    {
        mat_J[ij] = J_ij;
    }
}

auto
computeEriSSSS(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> CDenseMatrix
{
    const auto m_pi = mathconst::getPiValue();

    auto rank = mpi::rank(MPI_COMM_WORLD);

    if (rank % 8 == 0) cudaSafe(hipSetDevice(4));
    if (rank % 8 == 1) cudaSafe(hipSetDevice(5));
    if (rank % 8 == 2) cudaSafe(hipSetDevice(2));
    if (rank % 8 == 3) cudaSafe(hipSetDevice(3));
    if (rank % 8 == 4) cudaSafe(hipSetDevice(6));
    if (rank % 8 == 5) cudaSafe(hipSetDevice(7));
    if (rank % 8 == 6) cudaSafe(hipSetDevice(0));
    if (rank % 8 == 7) cudaSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Preparation");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("gpu::computeEriSSSS: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    // GTO info on device

    int64_t max_ncgtos = 0, max_npgtos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        max_ncgtos = std::max(ncgtos, max_ncgtos);
        max_npgtos = std::max(npgtos, max_npgtos);
    }

    double* d_gto_info;

    cudaSafe(hipMalloc(&d_gto_info, 5 * max_ncgtos * max_npgtos * sizeof(double)));

    // For now, only process one GTO block

    // TODO: take care of more GTO blocks

    auto gto_info = gpu::getGtoInfo(gto_blocks[0]);

    cudaSafe(hipMemcpy(d_gto_info, gto_info.data(), gto_info.size() * sizeof(double), hipMemcpyHostToDevice));

    auto ncgtos = gto_blocks[0].getNumberOfBasisFunctions();

    auto npgtos = gto_blocks[0].getNumberOfPrimitives();

    auto n_primitives = npgtos * ncgtos;

    std::vector<std::tuple<double, int64_t, int64_t, double>> mat_Q_D_sorted;

    auto dens_ptr = densityMatrix.alphaDensity(0);

    // For now, only process s-orbitals

    // TODO: go through spherical harmonics components

    for (int64_t i = 0; i < ncgtos * npgtos; i++)
    {
        const auto a_i = gto_info[i + n_primitives * 0];
        const auto c_i = gto_info[i + n_primitives * 1];
        const auto x_i = gto_info[i + n_primitives * 2];
        const auto y_i = gto_info[i + n_primitives * 3];
        const auto z_i = gto_info[i + n_primitives * 4];

        // const auto i_pgto = i / ncgtos;
        const auto i_cgto = i % ncgtos;

        for (int64_t j = 0; j <= i; j++)
        {
            const auto a_j = gto_info[j + n_primitives * 0];
            const auto c_j = gto_info[j + n_primitives * 1];
            const auto x_j = gto_info[j + n_primitives * 2];
            const auto y_j = gto_info[j + n_primitives * 3];
            const auto z_j = gto_info[j + n_primitives * 4];

            // const auto j_pgto = j / ncgtos;
            const auto j_cgto = j % ncgtos;

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij = c_i * c_j * std::pow(m_pi / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            // Note: t = 0.0 and F0_t = 1.0
            const auto eri_ijij = S_ij * S_ij * std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / m_pi);

            const auto sqrt_ijij = std::sqrt(eri_ijij);

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

            mat_Q_D_sorted.push_back(std::make_tuple(sqrt_ijij, i, j, D_ij));
        }
    }

    std::sort(mat_Q_D_sorted.begin(), mat_Q_D_sorted.end());
    std::reverse(mat_Q_D_sorted.begin(), mat_Q_D_sorted.end());

    // sorted Q, D, and indices on host

    auto n_prim_pairs = n_primitives * (n_primitives + 1) / 2;

    std::vector<double> mat_Q(n_prim_pairs);
    std::vector<double> mat_D(n_prim_pairs);
    std::vector<double> mat_J(n_prim_pairs);

    std::vector<uint32_t> first_inds(n_prim_pairs);
    std::vector<uint32_t> second_inds(n_prim_pairs);

    for (int64_t ij = 0; ij < n_prim_pairs; ij++)
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

    // sorted Q, D, and indices on device

    double *  d_mat_J, *d_mat_D, *d_mat_Q;
    uint32_t *d_first_inds, *d_second_inds;

    cudaSafe(hipMalloc(&d_mat_Q, n_prim_pairs * sizeof(double)));
    cudaSafe(hipMalloc(&d_mat_D, n_prim_pairs * sizeof(double)));
    cudaSafe(hipMalloc(&d_mat_J, n_prim_pairs * sizeof(double)));

    cudaSafe(hipMalloc(&d_first_inds, n_prim_pairs * sizeof(uint32_t)));
    cudaSafe(hipMalloc(&d_second_inds, n_prim_pairs * sizeof(uint32_t)));

    cudaSafe(hipMemcpy(d_mat_D, mat_D.data(), n_prim_pairs * sizeof(double), hipMemcpyHostToDevice));
    cudaSafe(hipMemcpy(d_mat_Q, mat_Q.data(), n_prim_pairs * sizeof(double), hipMemcpyHostToDevice));

    cudaSafe(hipMemcpy(d_first_inds, first_inds.data(), n_prim_pairs * sizeof(uint32_t), hipMemcpyHostToDevice));
    cudaSafe(hipMemcpy(d_second_inds, second_inds.data(), n_prim_pairs * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("Preparation");

    timer.start("J computation");

    // compute J

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((n_prim_pairs + threads_per_block.x - 1) / threads_per_block.x, 1);

    hipLaunchKernelGGL(cuda_mat_J_PI, num_blocks, threads_per_block, 0, 0, 
        d_mat_J, d_gto_info, static_cast<uint32_t>(n_primitives), d_mat_D, d_mat_Q, d_first_inds, d_second_inds, static_cast<uint32_t>(n_prim_pairs));

    cudaSafe(hipMemcpy(mat_J.data(), d_mat_J, n_prim_pairs * sizeof(double), hipMemcpyDeviceToHost));

    timer.stop("J computation");

    timer.start("J to Fock");

    // put J back to Fock

    CDenseMatrix mat_Fock(naos, naos);

    mat_Fock.zero();

    for (int64_t ij = 0; ij < n_prim_pairs; ij++)
    {
        const auto i = first_inds[ij];
        const auto j = second_inds[ij];

        const auto i_cgto = i % ncgtos;
        const auto j_cgto = j % ncgtos;

        mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij];

        if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij];
    }

    timer.stop("J to Fock");

    cudaSafe(hipFree(d_gto_info));
    cudaSafe(hipFree(d_mat_Q));
    cudaSafe(hipFree(d_mat_D));
    cudaSafe(hipFree(d_mat_J));
    cudaSafe(hipFree(d_first_inds));
    cudaSafe(hipFree(d_second_inds));

    timer.stop("Total timing");

    std::cout << "\nTiming of ERIs on GPU\n";
    std::cout << "-----------------------\n";
    std::cout << timer.getSummary() << std::endl;

    return mat_Fock;
}

}  // namespace gpu
