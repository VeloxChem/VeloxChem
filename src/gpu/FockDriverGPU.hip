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
#include <hipblas.h>
//#include <hipsolver.h>
#include <magma_v2.h>

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "ScreeningData.hpp"
#include "EriScreener.hpp"
#include "BoysFuncTable.hpp"
#include "DenseLinearAlgebra.hpp"
#include "FockDriverGPU.hpp"
#include "EriCoulombExchange.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MpiFunc.hpp"
#include "MultiTimer.hpp"
#include "OneElectronIntegrals.hpp"
#include "StringFormat.hpp"

#define TILE_DIM 16

#define MATH_CONST_PI 3.14159265358979323846

#define MATH_CONST_HALF_SQRT_PI 0.88622692545275794096

#define hipSafe(e)                                                                                                       \
    {                                                                                                                    \
        hipError_t err = (e);                                                                                            \
        if (err != hipSuccess)                                                                                           \
        {                                                                                                                \
            std::cerr << "HIP error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl;  \
            std::exit(EXIT_FAILURE);                                                                                     \
        }                                                                                                                \
    }

#define hipblasSafe(e)                                                                            \
    {                                                                                             \
        hipblasStatus_t err = (e);                                                                \
        if (err != HIPBLAS_STATUS_SUCCESS)                                                        \
        {                                                                                         \
            std::cerr << "hipBLAS error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                              \
        }                                                                                         \
    }

#define magmaSafe(e)                                                                                                   \
    {                                                                                                                  \
        magma_int_t err = (e);                                                                                         \
        if (err != MAGMA_SUCCESS) {                                                                                    \
            std::cerr << "MAGMA error in " << __FILE__ << ":" << __LINE__ << ": " << magma_strerror(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                   \
        }                                                                                                              \
    }

/*
#define hipsolverSafe(e)                                                                            \
    {                                                                                               \
        hipsolverStatus_t err = (e);                                                                \
        if (err != HIPSOLVER_STATUS_SUCCESS)                                                        \
        {                                                                                           \
            std::cerr << "hipSolver error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                                \
        }                                                                                           \
    }
*/

namespace gpu {  // gpu namespace

auto
computeOneElectronIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::vector<CDenseMatrix> STV_matrices(3);

    STV_matrices[0] = CDenseMatrix(naos, naos);
    STV_matrices[1] = CDenseMatrix(naos, naos);
    STV_matrices[2] = CDenseMatrix(naos, naos);

    for (int64_t ind = 0; ind < naos * naos; ind++)
    {
        STV_matrices[0].values()[ind] = 0.0;
        STV_matrices[1].values()[ind] = 0.0;
        STV_matrices[2].values()[ind] = 0.0;
    }

    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    std::vector<CDenseMatrix> S_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> T_matrices(num_gpus_per_node);
    std::vector<CDenseMatrix> V_matrices(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        S_matrices[gpu_id] = CDenseMatrix(naos, naos);
        T_matrices[gpu_id] = CDenseMatrix(naos, naos);
        V_matrices[gpu_id] = CDenseMatrix(naos, naos);
    }

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    // std::stringstream ss;
    // ss << "rank " << rank << ",  thread " << thread_id << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    if (gpu_rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (gpu_rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (gpu_rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (gpu_rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (gpu_rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (gpu_rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (gpu_rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (gpu_rank % 8 == 7) hipSafe(hipSetDevice(1));

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    hipSafe(hipMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;
    uint32_t* d_s_prim_aoinds;

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_s_prim_aoinds, s_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;
    uint32_t* d_p_prim_aoinds;

    hipSafe(hipMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_p_prim_aoinds, p_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // GTO block pairs

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, pp_prim_pair_count_local});

    std::vector<double> mat_S(max_prim_pair_count_local);
    std::vector<double> mat_T(max_prim_pair_count_local);
    std::vector<double> mat_V(max_prim_pair_count_local);

    // S on device

    double *d_mat_S, *d_mat_T, *d_mat_V;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;

    hipSafe(hipMalloc(&d_mat_S, max_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_T, max_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_V, max_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    for (int64_t ind = 0; ind < naos * naos; ind++)
    {
        S_matrices[gpu_id].values()[ind] = 0.0;
        T_matrices[gpu_id].values()[ind] = 0.0;
        V_matrices[gpu_id].values()[ind] = 0.0;
    }

    auto& mat_overlap = S_matrices[gpu_id];
    auto& mat_kinetic_energy = T_matrices[gpu_id];
    auto& mat_nuclear_potential = V_matrices[gpu_id];

    const auto mol_charges = molecule.getCharges();
    const auto mol_coords = molecule.getCoordinates(std::string("BOHR"));
    const auto natoms = molecule.getNumberOfAtoms();

    std::vector<double> points_info(natoms * 4);

    for (int64_t a = 0; a < natoms; a++)
    {
        points_info[a + natoms * 0] = mol_coords[a][0];
        points_info[a + natoms * 1] = mol_coords[a][1];
        points_info[a + natoms * 2] = mol_coords[a][2];
        points_info[a + natoms * 3] = mol_charges[a];
    }

    double *d_points_info;

    hipSafe(hipMalloc(&d_points_info, natoms * 4 * sizeof(double)));

    hipSafe(hipMemcpy(d_points_info, points_info.data(), natoms * 4 * sizeof(double), hipMemcpyHostToDevice));

    // compute S

    // S: SS

    if (ss_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergySS, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialSS, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_overlap.row(i_cgto)[j_cgto] += mat_S[ij];

            if (i != j) mat_overlap.row(j_cgto)[i_cgto] += mat_S[ij];

            mat_kinetic_energy.row(i_cgto)[j_cgto] += mat_T[ij];

            if (i != j) mat_kinetic_energy.row(j_cgto)[i_cgto] += mat_T[ij];

            mat_nuclear_potential.row(i_cgto)[j_cgto] += mat_V[ij];

            if (i != j) mat_nuclear_potential.row(j_cgto)[i_cgto] += mat_V[ij];
        }
    }

    // S: SP

    if (sp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergySP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialSP, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_overlap.row(i_cgto)[j_cgto] += mat_S[ij];
            mat_overlap.row(j_cgto)[i_cgto] += mat_S[ij];

            mat_kinetic_energy.row(i_cgto)[j_cgto] += mat_T[ij];
            mat_kinetic_energy.row(j_cgto)[i_cgto] += mat_T[ij];

            mat_nuclear_potential.row(i_cgto)[j_cgto] += mat_V[ij];
            mat_nuclear_potential.row(j_cgto)[i_cgto] += mat_V[ij];
        }
    }

    // S: PP

    if (pp_prim_pair_count_local > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::computeOverlapAndKineticEnergyPP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local));

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        hipLaunchKernelGGL(gpu::computeNuclearPotentialPP, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_V.data(), d_mat_V, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_overlap.row(i_cgto)[j_cgto] += mat_S[ij];

            if (i != j) mat_overlap.row(j_cgto)[i_cgto] += mat_S[ij];

            mat_kinetic_energy.row(i_cgto)[j_cgto] += mat_T[ij];

            if (i != j) mat_kinetic_energy.row(j_cgto)[i_cgto] += mat_T[ij];

            mat_nuclear_potential.row(i_cgto)[j_cgto] += mat_V[ij];

            if (i != j) mat_nuclear_potential.row(j_cgto)[i_cgto] += mat_V[ij];
        }
    }

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_s_prim_aoinds));

    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_p_prim_aoinds));

    hipSafe(hipFree(d_points_info));

    hipSafe(hipFree(d_mat_S));
    hipSafe(hipFree(d_mat_T));
    hipSafe(hipFree(d_mat_V));

    hipSafe(hipFree(d_ss_first_inds_local));
    hipSafe(hipFree(d_ss_second_inds_local));
    hipSafe(hipFree(d_sp_first_inds_local));
    hipSafe(hipFree(d_sp_second_inds_local));
    hipSafe(hipFree(d_pp_first_inds_local));
    hipSafe(hipFree(d_pp_second_inds_local));

    }}

    auto p_mat_S = STV_matrices[0].values();
    auto p_mat_T = STV_matrices[1].values();
    auto p_mat_V = STV_matrices[2].values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_overlap = S_matrices[gpu_id].values();
        auto p_mat_kinetic_energy = T_matrices[gpu_id].values();
        auto p_mat_nuclear_potential = V_matrices[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_S[ind] += p_mat_overlap[ind];
            p_mat_T[ind] += p_mat_kinetic_energy[ind];
            p_mat_V[ind] += p_mat_nuclear_potential[ind];
        }
    }

    return STV_matrices;
}

auto
computeFockOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix, CScreeningData& screening) -> CDenseMatrix
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    auto nthreads = omp_get_max_threads();
    auto num_gpus_per_node = screening.getNumGpusPerNode();
    auto num_threads_per_gpu = nthreads / num_gpus_per_node;

    auto gpu_rank = rank * num_gpus_per_node;

    // std::stringstream ss;
    // ss << "\n(Prep.) rank " << rank << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    if (gpu_rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (gpu_rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (gpu_rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (gpu_rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (gpu_rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (gpu_rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (gpu_rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (gpu_rank % 8 == 7) hipSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Prep. blocks");

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);
    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
    }

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    timer.stop("Prep. blocks");

    timer.start("Prep. sortQD");

    auto dens_ptr = densityMatrix.alphaDensity(0);

    screening.sortQD(s_prim_count, p_prim_count, s_prim_aoinds, p_prim_aoinds, naos, dens_ptr);

    timer.stop("Prep. sortQD");

    timer.start("Prep. Q_prime");

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    // TODO distribute computation of Q_prime

    double *d_matrix_A, *d_matrix_B, *d_matrix_C;

    auto mat_full = screening.get_mat_Q_full(s_prim_count, p_prim_count);

    hipSafe(hipMalloc(&d_matrix_A, mat_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_B, mat_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_C, mat_full.getNumberOfElements() * sizeof(double)));

    hipSafe(hipMemcpy(d_matrix_A, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));

    mat_full = screening.get_mat_D_abs_full(s_prim_count, p_prim_count, s_prim_aoinds, p_prim_aoinds, naos, dens_ptr);

    hipSafe(hipMemcpy(d_matrix_B, mat_full.values(), mat_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));

    const auto all_prim_count = s_prim_count + p_prim_count * 3;

    // TODO double check the criteria
    if (all_prim_count < 8192)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((all_prim_count + threads_per_block.x - 1) / threads_per_block.x, (all_prim_count + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::matmulAB, num_blocks, threads_per_block, 0, 0, d_matrix_C, d_matrix_A, d_matrix_B,
                static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count));

        hipSafe(hipDeviceSynchronize());

        hipLaunchKernelGGL(gpu::matmulAB, num_blocks, threads_per_block, 0, 0, d_matrix_B, d_matrix_C, d_matrix_A,
                static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count));
    }
    else
    {
        hipblasHandle_t handle;
        hipblasSafe(hipblasCreate(&handle));

        double alpha = 1.0, beta = 0.0;

        auto n = static_cast<int32_t>(all_prim_count);

        // compute A^T * (B^T * A^T) since hipblas is column-major
        hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, n, n, n, &alpha, d_matrix_B, n, d_matrix_A, n, &beta, d_matrix_C, n));

        hipSafe(hipDeviceSynchronize());

        hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, n, n, n, &alpha, d_matrix_A, n, d_matrix_C, n, &beta, d_matrix_B, n));

        hipblasSafe(hipblasDestroy(handle));
    }

    hipSafe(hipMemcpy(mat_full.values(), d_matrix_B, mat_full.getNumberOfElements() * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_matrix_A));
    hipSafe(hipFree(d_matrix_B));
    hipSafe(hipFree(d_matrix_C));

    timer.stop("Prep. Q_prime");

    timer.start("Prep. preLinK");

    // auto& Q_prime = mat_full;
    const double Q_prime_thresh = 1.0e-7;

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, mat_full, Q_prime_thresh);

    timer.stop("Prep. preLinK");

    std::string errnaos("gpu::computeCoulombFock: Inconsistent number of AOs");
    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    CDenseMatrix mat_Fock_sum (naos, naos);

    for (int64_t ind = 0; ind < naos * naos; ind++)
    {
        mat_Fock_sum.values()[ind] = 0.0;
    }

    std::vector<CDenseMatrix> mat_Fock_omp(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        mat_Fock_omp[gpu_id] = CDenseMatrix(naos, naos);
    }

    std::cout << "\nTiming of prep. on rank " << rank << "\n";
    std::cout << "-------------------------\n";
    std::cout << timer.getSummary() << "\n";

#pragma omp parallel
    {
    auto thread_id = omp_get_thread_num();

    if (thread_id % num_threads_per_gpu == 0)
    {
    auto gpu_id = thread_id / num_threads_per_gpu;
    auto gpu_rank = gpu_id + rank * num_gpus_per_node;
    auto gpu_count = nnodes * num_gpus_per_node;

    // std::stringstream ss;
    // ss << "rank " << rank << ",  thread " << thread_id << ",  gpu_rank " << gpu_rank << "\n";
    // std::cout << ss.str();

    if (gpu_rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (gpu_rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (gpu_rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (gpu_rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (gpu_rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (gpu_rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (gpu_rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (gpu_rank % 8 == 7) hipSafe(hipSetDevice(1));

    CMultiTimer timer;

    timer.start("Total timing");

    timer.start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    double* d_boys_func_ft;

    hipSafe(hipMalloc(&d_boys_func_ft, boys_func_ft.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));

    timer.stop("Boys func. prep.");

    timer.start("GTO block prep.");

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
    }

    // S gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    double*   d_s_prim_info;
    uint32_t* d_s_prim_aoinds;

    hipSafe(hipMalloc(&d_s_prim_info, s_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_s_prim_aoinds, s_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_s_prim_aoinds, s_prim_aoinds.data(), s_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    // P gto block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    double*   d_p_prim_info;
    uint32_t* d_p_prim_aoinds;

    hipSafe(hipMalloc(&d_p_prim_info, p_prim_info.size() * sizeof(double)));
    hipSafe(hipMalloc(&d_p_prim_aoinds, p_prim_aoinds.size() * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_p_prim_aoinds, p_prim_aoinds.data(), p_prim_aoinds.size() * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("GTO block prep.");

    timer.start("GTO block pair prep.");

    // GTO block pairs

    auto dens_ptr = densityMatrix.alphaDensity(0);

    // TODO distribute computation of Q matrices

    // timer.start("  D matrices");

    // screening.sortQD(s_prim_count, p_prim_count, s_prim_aoinds, p_prim_aoinds, naos, dens_ptr);

    // timer.stop("  D matrices");

    timer.stop("GTO block pair prep.");

    timer.start("Coulomb prep.");

    const auto ss_mat_Q_orig = screening.getQMatrixSS();
    const auto sp_mat_Q_orig = screening.getQMatrixSP();
    const auto pp_mat_Q_orig = screening.getQMatrixPP();

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(gpu_id);
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(gpu_id);
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(gpu_id);

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(gpu_id);
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(gpu_id);
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(gpu_id);

    const auto& ss_mat_Q_local = screening.get_ss_mat_Q_local(gpu_id);
    const auto& sp_mat_Q_local = screening.get_sp_mat_Q_local(gpu_id);
    const auto& pp_mat_Q_local = screening.get_pp_mat_Q_local(gpu_id);

    const auto& ss_first_inds = screening.get_ss_first_inds();
    const auto& sp_first_inds = screening.get_sp_first_inds();
    const auto& pp_first_inds = screening.get_pp_first_inds();

    const auto& ss_second_inds = screening.get_ss_second_inds();
    const auto& sp_second_inds = screening.get_sp_second_inds();
    const auto& pp_second_inds = screening.get_pp_second_inds();

    const auto& ss_mat_Q = screening.get_ss_mat_Q(); 
    const auto& sp_mat_Q = screening.get_sp_mat_Q(); 
    const auto& pp_mat_Q = screening.get_pp_mat_Q(); 

    const auto& ss_mat_D = screening.get_ss_mat_D(); 
    const auto& sp_mat_D = screening.get_sp_mat_D(); 
    const auto& pp_mat_D = screening.get_pp_mat_D(); 

    const auto ss_max_D = screening.get_ss_max_D(); 
    const auto sp_max_D = screening.get_sp_max_D(); 
    const auto pp_max_D = screening.get_pp_max_D(); 

    const auto ss_prim_pair_count = static_cast<int64_t>(ss_first_inds.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sp_first_inds.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(pp_first_inds.size());

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, pp_prim_pair_count});

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());

    const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local, pp_prim_pair_count_local});

    std::vector<double> mat_J(max_prim_pair_count_local);

    // sorted Q, D, and indices on device

    double *d_mat_J, *d_mat_D, *d_ss_mat_Q, *d_sp_mat_Q, *d_pp_mat_Q;

    uint32_t *d_ss_first_inds, *d_ss_second_inds;
    uint32_t *d_sp_first_inds, *d_sp_second_inds;
    uint32_t *d_pp_first_inds, *d_pp_second_inds;

    double *d_ss_mat_Q_local, *d_sp_mat_Q_local, *d_pp_mat_Q_local;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;

    hipSafe(hipMalloc(&d_mat_D, max_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_J, max_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_mat_Q, ss_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_mat_Q, sp_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_mat_Q, pp_prim_pair_count * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds, ss_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds, ss_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds, sp_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds, sp_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds, pp_prim_pair_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds, pp_prim_pair_count * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_ss_mat_Q_local, ss_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_sp_mat_Q_local, sp_prim_pair_count_local * sizeof(double)));
    hipSafe(hipMalloc(&d_pp_mat_Q_local, pp_prim_pair_count_local * sizeof(double)));

    hipSafe(hipMalloc(&d_ss_first_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_ss_second_inds_local, ss_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_sp_first_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_sp_second_inds_local, sp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pp_first_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pp_second_inds_local, pp_prim_pair_count_local * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_ss_mat_Q, ss_mat_Q.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_mat_Q, sp_mat_Q.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_mat_Q, pp_mat_Q.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_first_inds, ss_first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds, sp_first_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds, pp_first_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_mat_Q_local, ss_mat_Q_local.data(), ss_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_mat_Q_local, sp_mat_Q_local.data(), sp_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_mat_Q_local, pp_mat_Q_local.data(), pp_prim_pair_count_local * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_first_inds_local, ss_first_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds_local, ss_second_inds_local.data(), ss_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds_local, sp_first_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds_local, sp_second_inds_local.data(), sp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds_local, pp_first_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds_local, pp_second_inds_local.data(), pp_prim_pair_count_local * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("Coulomb prep.");

    for (int64_t ind = 0; ind < naos * naos; ind++)
    {
        mat_Fock_omp[gpu_id].values()[ind] = 0.0;
    }

    // auto& mat_Fock = mat_Fock_omp[gpu_id];

    timer.start("J computation");

    // compute J

    // J: S-S block

    if (ss_prim_pair_count_local > 0)
    {
        timer.start("  J block SS");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(ss_prim_pair_count_local));

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SS|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSSS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_mat_D,
                                                     d_ss_mat_Q_local,
                                                     d_ss_mat_Q,
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        // J: (SS|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSSP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_ss_mat_Q_local,
                                                     d_sp_mat_Q,
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        // J: (SS|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSSPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_ss_mat_Q_local,
                                                     d_pp_mat_Q,
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }

        timer.stop("  J block SS");
    }

    // J: S-P block

    if (sp_prim_pair_count_local > 0)
    {
        timer.start("  J block SP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(sp_prim_pair_count_local));

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (SP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPSS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_sp_mat_Q_local,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        }

        // J: (SP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPSP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_sp_mat_Q_local,
                                                     d_sp_mat_Q,
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local),
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        // J: (SP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockSPPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_sp_mat_Q_local,
                                                     d_pp_mat_Q,
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;
            mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }

        timer.stop("  J block SP");
    }

    // J: P-P block

    if (pp_prim_pair_count_local > 0)
    {
        timer.start("  J block PP");

        // zeroize J on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_J, static_cast<uint32_t>(pp_prim_pair_count_local));

        // set up thread blocks for J

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        // J: (PP|SS)
        //     **

        if (ss_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPSS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_pp_mat_Q_local,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        // J: (PP|SP)
        //     **

        if (sp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPSP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_sp_mat_Q,
                                                     d_pp_mat_Q_local,
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        // J: (PP|PP)
        //     **

        if (pp_prim_pair_count > 0)
        {
            hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

            hipLaunchKernelGGL(gpu::computeCoulombFockPPPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_J,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_mat_D,
                                                     d_pp_mat_Q_local,
                                                     d_pp_mat_Q,
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);
        }

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock_omp[gpu_id].row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock_omp[gpu_id].row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }

        timer.stop("  J block PP");
    }

    hipSafe(hipFree(d_mat_D));
    hipSafe(hipFree(d_mat_J));

    hipSafe(hipFree(d_ss_mat_Q));
    hipSafe(hipFree(d_sp_mat_Q));
    hipSafe(hipFree(d_pp_mat_Q));

    hipSafe(hipFree(d_ss_first_inds));
    hipSafe(hipFree(d_ss_second_inds));
    hipSafe(hipFree(d_sp_first_inds));
    hipSafe(hipFree(d_sp_second_inds));
    hipSafe(hipFree(d_pp_first_inds));
    hipSafe(hipFree(d_pp_second_inds));

    hipSafe(hipFree(d_ss_mat_Q_local));
    hipSafe(hipFree(d_sp_mat_Q_local));
    hipSafe(hipFree(d_pp_mat_Q_local));

    hipSafe(hipFree(d_ss_first_inds_local));
    hipSafe(hipFree(d_ss_second_inds_local));
    hipSafe(hipFree(d_sp_first_inds_local));
    hipSafe(hipFree(d_sp_second_inds_local));
    hipSafe(hipFree(d_pp_first_inds_local));
    hipSafe(hipFree(d_pp_second_inds_local));

    timer.stop("J computation");

    timer.start("Exchange prep.");

    // K preparation

    const auto& mat_Q_for_K_ss = screening.get_mat_Q_for_K_ss();
    const auto& mat_Q_for_K_sp = screening.get_mat_Q_for_K_sp();
    const auto& mat_Q_for_K_ps = screening.get_mat_Q_for_K_ps();
    const auto& mat_Q_for_K_pp = screening.get_mat_Q_for_K_pp();

    const auto& density_inds_for_K_ss = screening.get_density_inds_for_K_ss();
    const auto& density_inds_for_K_sp = screening.get_density_inds_for_K_sp();
    const auto& density_inds_for_K_ps = screening.get_density_inds_for_K_ps();
    const auto& density_inds_for_K_pp = screening.get_density_inds_for_K_pp();

    const auto& pair_inds_i_for_K_ss = screening.get_local_pair_inds_i_for_K_ss(gpu_id);
    const auto& pair_inds_k_for_K_ss = screening.get_local_pair_inds_k_for_K_ss(gpu_id);

    const auto& pair_inds_i_for_K_sp = screening.get_local_pair_inds_i_for_K_sp(gpu_id);
    const auto& pair_inds_k_for_K_sp = screening.get_local_pair_inds_k_for_K_sp(gpu_id);

    const auto& pair_inds_i_for_K_pp = screening.get_local_pair_inds_i_for_K_pp(gpu_id);
    const auto& pair_inds_k_for_K_pp = screening.get_local_pair_inds_k_for_K_pp(gpu_id);

    auto pair_inds_count_for_K_ss = static_cast<uint32_t>(pair_inds_i_for_K_ss.size());
    auto pair_inds_count_for_K_sp = static_cast<uint32_t>(pair_inds_i_for_K_sp.size());
    auto pair_inds_count_for_K_pp = static_cast<uint32_t>(pair_inds_i_for_K_pp.size());

    double*   d_mat_K;
    uint32_t *d_pair_inds_i_for_K_ss, *d_pair_inds_k_for_K_ss;
    uint32_t *d_pair_inds_i_for_K_sp, *d_pair_inds_k_for_K_sp;
    uint32_t *d_pair_inds_i_for_K_pp, *d_pair_inds_k_for_K_pp;

    const auto max_pair_inds_count = std::max({pair_inds_count_for_K_ss, pair_inds_count_for_K_sp, pair_inds_count_for_K_pp});

    std::vector<double> mat_K(max_pair_inds_count);

    hipSafe(hipMalloc(&d_mat_K, max_pair_inds_count * sizeof(double)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_pp, pair_inds_count_for_K_pp * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_pp, pair_inds_i_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_pp, pair_inds_k_for_K_pp.data(), pair_inds_count_for_K_pp * sizeof(uint32_t), hipMemcpyHostToDevice));

    double* d_mat_D_full_AO;
    double *d_mat_Q_for_K_ss, *d_mat_Q_for_K_sp, *d_mat_Q_for_K_ps, *d_mat_Q_for_K_pp;
    uint32_t *d_density_inds_for_K_ss, *d_density_inds_for_K_sp, *d_density_inds_for_K_ps, *d_density_inds_for_K_pp;

    hipSafe(hipMalloc(&d_mat_D_full_AO, naos * naos * sizeof(double)));

    hipSafe(hipMalloc(&d_mat_Q_for_K_ss, s_prim_count * s_prim_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Q_for_K_sp, s_prim_count * p_prim_count * 3 * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Q_for_K_ps, p_prim_count * 3 * s_prim_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Q_for_K_pp, p_prim_count * 3 * p_prim_count * 3 * sizeof(double)));

    hipSafe(hipMalloc(&d_density_inds_for_K_ss, s_prim_count * s_prim_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_density_inds_for_K_sp, s_prim_count * p_prim_count * 3 * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_density_inds_for_K_ps, p_prim_count * 3 * s_prim_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_density_inds_for_K_pp, p_prim_count * 3 * p_prim_count * 3 * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_mat_D_full_AO, dens_ptr, naos * naos * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_mat_Q_for_K_ss, mat_Q_for_K_ss.data(), s_prim_count * s_prim_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q_for_K_sp, mat_Q_for_K_sp.data(), s_prim_count * p_prim_count * 3 * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q_for_K_ps, mat_Q_for_K_ps.data(), p_prim_count * 3 * s_prim_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q_for_K_pp, mat_Q_for_K_pp.data(), p_prim_count * 3 * p_prim_count * 3 * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_density_inds_for_K_ss, density_inds_for_K_ss.data(), s_prim_count * s_prim_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_density_inds_for_K_sp, density_inds_for_K_sp.data(), s_prim_count * p_prim_count * 3 * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_density_inds_for_K_ps, density_inds_for_K_ps.data(), p_prim_count * 3 * s_prim_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_density_inds_for_K_pp, density_inds_for_K_pp.data(), p_prim_count * 3 * p_prim_count * 3 * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("Exchange prep.");

    timer.start("K computation");

    // compute K

    // K: S-S block

    if (pair_inds_count_for_K_ss > 0)
    {
        timer.start("  K block SS");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_ss + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_ss));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3(pair_inds_count_for_K_ss, 1);

        // K: (SS|SS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSSS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_ss,
                                                     d_pair_inds_k_for_K_ss,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_ss),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     ss_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ss,
                                                     d_density_inds_for_K_ss,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|SP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSSP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_ss,
                                                     d_pair_inds_k_for_K_ss,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_ss),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ss,
                                                     d_mat_Q_for_K_sp,
                                                     d_density_inds_for_K_ss,
                                                     d_density_inds_for_K_sp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|SS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPSS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_ss,
                                                     d_pair_inds_k_for_K_ss,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_ss),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ss,
                                                     d_mat_Q_for_K_sp,
                                                     d_density_inds_for_K_ss,
                                                     d_density_inds_for_K_sp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|SP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPSP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_ss,
                                                     d_pair_inds_k_for_K_ss,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_ss),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_sp,
                                                     d_density_inds_for_K_sp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock_omp[gpu_id].row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }

        timer.stop("  K block SS");
    }

    // K: S-P block

    if (pair_inds_count_for_K_sp > 0)
    {
        timer.start("  K block SP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_sp + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_sp));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3(pair_inds_count_for_K_sp, 1);

        // K: (SS|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSPS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_sp,
                                                     d_pair_inds_k_for_K_sp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_sp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     ss_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ss,
                                                     d_mat_Q_for_K_ps,
                                                     d_density_inds_for_K_ss,
                                                     d_density_inds_for_K_ps,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SS|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSSPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_sp,
                                                     d_pair_inds_k_for_K_sp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_sp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ss,
                                                     d_mat_Q_for_K_pp,
                                                     d_density_inds_for_K_ss,
                                                     d_density_inds_for_K_pp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPPS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_sp,
                                                     d_pair_inds_k_for_K_sp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_sp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_sp,
                                                     d_mat_Q_for_K_ps,
                                                     d_density_inds_for_K_sp,
                                                     d_density_inds_for_K_ps,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (SP|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockSPPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_sp,
                                                     d_pair_inds_k_for_K_sp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_sp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_sp,
                                                     d_mat_Q_for_K_pp,
                                                     d_density_inds_for_K_sp,
                                                     d_density_inds_for_K_pp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);
            mat_Fock_omp[gpu_id].row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }

        timer.stop("  K block SP");
    }

    // K: P-P block

    if (pair_inds_count_for_K_pp > 0)
    {
        timer.start("  K block PP");

        // zeroize K on device

        dim3 threads_per_block(TILE_DIM * TILE_DIM);

        dim3 num_blocks((pair_inds_count_for_K_pp + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(gpu::zeroData, num_blocks, threads_per_block, 0, 0, d_mat_K, static_cast<uint32_t>(pair_inds_count_for_K_pp));

        hipSafe(hipDeviceSynchronize());

        // set up thread blocks for K

        threads_per_block = dim3(TILE_DIM, TILE_DIM);

        num_blocks = dim3(pair_inds_count_for_K_pp, 1);

        // K: (PS|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockPSPS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_pp,
                                                     d_pair_inds_k_for_K_pp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_pp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     ss_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ps,
                                                     d_density_inds_for_K_ps,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (PS|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockPSPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_pp,
                                                     d_pair_inds_k_for_K_pp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_pp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ps,
                                                     d_mat_Q_for_K_pp,
                                                     d_density_inds_for_K_ps,
                                                     d_density_inds_for_K_pp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (PP|PS)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPS, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_pp,
                                                     d_pair_inds_k_for_K_pp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_pp),
                                                     d_s_prim_info,
                                                     d_s_prim_aoinds,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_ps,
                                                     d_mat_Q_for_K_pp,
                                                     d_density_inds_for_K_ps,
                                                     d_density_inds_for_K_pp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipDeviceSynchronize());

        // K: (PP|PP)
        //     *  *

        hipLaunchKernelGGL(gpu::computeExchangeFockPPPP, num_blocks, threads_per_block, 0, 0,
                                                     d_mat_K,
                                                     d_pair_inds_i_for_K_pp,
                                                     d_pair_inds_k_for_K_pp,
                                                     static_cast<uint32_t>(pair_inds_count_for_K_pp),
                                                     d_p_prim_info,
                                                     d_p_prim_aoinds,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D_full_AO,
                                                     d_mat_Q_for_K_pp,
                                                     d_density_inds_for_K_pp,
                                                     static_cast<uint32_t>(naos),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock_omp[gpu_id].row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock_omp[gpu_id].row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }

        timer.stop("  K block PP");
    }

    timer.stop("K computation");

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_s_prim_aoinds));

    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_p_prim_aoinds));

    hipSafe(hipFree(d_mat_K));

    hipSafe(hipFree(d_pair_inds_i_for_K_ss));
    hipSafe(hipFree(d_pair_inds_k_for_K_ss));
    hipSafe(hipFree(d_pair_inds_i_for_K_sp));
    hipSafe(hipFree(d_pair_inds_k_for_K_sp));
    hipSafe(hipFree(d_pair_inds_i_for_K_pp));
    hipSafe(hipFree(d_pair_inds_k_for_K_pp));

    hipSafe(hipFree(d_mat_D_full_AO));

    hipSafe(hipFree(d_mat_Q_for_K_ss));
    hipSafe(hipFree(d_mat_Q_for_K_sp));
    hipSafe(hipFree(d_mat_Q_for_K_ps));
    hipSafe(hipFree(d_mat_Q_for_K_pp));

    hipSafe(hipFree(d_density_inds_for_K_ss));
    hipSafe(hipFree(d_density_inds_for_K_sp));
    hipSafe(hipFree(d_density_inds_for_K_ps));
    hipSafe(hipFree(d_density_inds_for_K_pp));

    timer.stop("Total timing");

    std::cout << "\nTiming of ERIs on GPU " << gpu_rank << "\n";
    std::cout << "-----------------------\n";
    std::cout << timer.getSummary() << std::endl;

    }}

    auto p_mat_F = mat_Fock_sum.values();

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        auto p_mat_fock = mat_Fock_omp[gpu_id].values();

        for (int64_t ind = 0; ind < naos * naos; ind++)
        {
            p_mat_F[ind] += p_mat_fock[ind];
        }
    }

    return mat_Fock_sum;
}

auto
computeMatrixMultiplication(double* C, const double* A, const double* B, const int64_t nrows_A, const int64_t ncols_A, const int64_t ncols_B) -> void
{
    // TODO allow computeMatrixMultiplication on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::computeMatrixMultiplication: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(4));

    const auto nrows_B = ncols_A;

    double *d_A, *d_B, *d_C;

    hipSafe(hipMalloc(&d_A, nrows_A * ncols_A * sizeof(double)));
    hipSafe(hipMalloc(&d_B, nrows_B * ncols_B * sizeof(double)));
    hipSafe(hipMalloc(&d_C, nrows_A * ncols_B * sizeof(double)));

    hipSafe(hipMemcpy(d_A, A, nrows_A * ncols_A * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_B, B, nrows_B * ncols_B * sizeof(double), hipMemcpyHostToDevice));

    // TODO double check the criteria
    if (std::max({nrows_A, ncols_A, ncols_B}) < 8192)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ncols_B + threads_per_block.x - 1) / threads_per_block.x, (nrows_A + threads_per_block.y - 1) / threads_per_block.y);

        hipLaunchKernelGGL(gpu::matmulAB, num_blocks, threads_per_block, 0, 0, d_C, d_A, d_B,
                static_cast<uint32_t>(nrows_A), static_cast<uint32_t>(ncols_A), static_cast<uint32_t>(ncols_B));
    }
    else
    {
        hipblasHandle_t handle;
        hipblasSafe(hipblasCreate(&handle));

        double alpha = 1.0, beta = 0.0;

        auto m = static_cast<int32_t>(nrows_A);
        auto k = static_cast<int32_t>(ncols_A);
        auto n = static_cast<int32_t>(ncols_B);

        // compute C^T = B^T * A^T since hipblas is column-major
        hipblasSafe(hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, n, m, k, &alpha, d_B, n, d_A, k, &beta, d_C, n));

        hipblasSafe(hipblasDestroy(handle));
    }

    hipSafe(hipMemcpy(C, d_C, nrows_A * ncols_B * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_B));
    hipSafe(hipFree(d_C));
}

auto
diagonalizeMatrix(double* A, double* W, const int64_t nrows_A) -> void
{
    // TODO allow diagonalizeMatrix on non-master MPI rank

    auto rank = mpi::rank(MPI_COMM_WORLD);

    std::string errmpirank("gpu::diagonalizeMatrix: Should only be called from MPI master rank");
    errors::assertMsgCritical(rank == mpi::master(), errmpirank);

    hipSafe(hipSetDevice(4));

    magmaSafe(magma_init());

    magma_setdevice(4);

    magma_int_t n = static_cast<magma_int_t>(nrows_A);
    //magma_int_t ldda = magma_roundup(n, 32);

    double *d_A;
    //hipSafe(hipMalloc(&d_A, n * ldda * sizeof(double)));
    //hipblasSafe(hipblasSetMatrix(n, n, sizeof(double), A, n, d_A, ldda));
    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    auto nb = magma_get_dsytrd_nb(n);
    auto lwork = static_cast<magma_int_t>(std::max(2*n + n*nb, 1 + 6*n + 2*n*n));
    auto liwork = 3 + 5*n;

    double *wA, *work;
    magma_int_t *iwork;

    hipSafe(hipHostMalloc(&wA, n * n * sizeof(double)));
    hipSafe(hipHostMalloc(&work, lwork * sizeof(double)));
    hipSafe(hipHostMalloc(&iwork, liwork * sizeof(magma_int_t)));

    magma_int_t info;
    //magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, ldda, W, wA, n, work, lwork, iwork, liwork, &info);
    magma_dsyevd_gpu(MagmaVec, MagmaUpper, n, d_A, n, W, wA, n, work, lwork, iwork, liwork, &info);

    if (info != 0)
    {
        std::stringstream ss;
        ss << "gpu::diagonalizeMatrix: (magma error) " << magma_strerror(info);
        errors::assertMsgCritical(false, ss.str());
    }

    //hipblasSafe(hipblasGetMatrix(n, n, sizeof(double), d_A, ldda, A, n));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));

    hipSafe(hipHostFree(wA));
    hipSafe(hipHostFree(work));
    hipSafe(hipHostFree(iwork));

    magmaSafe(magma_finalize());
}

/*
auto
diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void
{
    hipSafe(hipSetDevice(4));

    auto n = static_cast<int32_t>(nrows_A);
    int32_t lwork;

    double *d_A, *d_D, *d_work;
    int32_t *d_info;

    hipSafe(hipMalloc(&d_A, n * n * sizeof(double)));
    hipSafe(hipMalloc(&d_D, n * sizeof(double)));
    hipSafe(hipMalloc(&d_info, sizeof(int32_t)));

    hipSafe(hipMemcpy(d_A, A, n * n * sizeof(double), hipMemcpyHostToDevice));

    hipsolverHandle_t handle;
    hipsolverSafe(hipsolverCreate(&handle));

    hipsolverSafe(hipsolverDsyevd_bufferSize(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, &lwork));

    hipSafe(hipMalloc(&d_work, lwork));

    hipsolverSafe(hipsolverDsyevd(handle, HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, n, d_A, n, d_D, d_work, lwork, d_info));

    hipsolverSafe(hipsolverDestroy(handle));

    // TODO check h_info
    int32_t h_info;

    hipSafe(hipMemcpy(&h_info, d_info, sizeof(int32_t), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(A, d_A, n * n * sizeof(double), hipMemcpyDeviceToHost));
    hipSafe(hipMemcpy(D, d_D, n * sizeof(double), hipMemcpyDeviceToHost));

    hipSafe(hipFree(d_A));
    hipSafe(hipFree(d_D));
    hipSafe(hipFree(d_info));
    hipSafe(hipFree(d_work));
}
*/

}  // namespace gpu
