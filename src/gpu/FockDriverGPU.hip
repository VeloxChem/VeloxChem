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

auto
computeOneElectronIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    if (rank % 8 == 0) hipSafe(hipSetDevice(4));
    if (rank % 8 == 1) hipSafe(hipSetDevice(5));
    if (rank % 8 == 2) hipSafe(hipSetDevice(2));
    if (rank % 8 == 3) hipSafe(hipSetDevice(3));
    if (rank % 8 == 4) hipSafe(hipSetDevice(6));
    if (rank % 8 == 5) hipSafe(hipSetDevice(7));
    if (rank % 8 == 6) hipSafe(hipSetDevice(0));
    if (rank % 8 == 7) hipSafe(hipSetDevice(1));

    std::vector<CDenseMatrix> STV_matrices(3);

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

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(); 
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(); 
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(); 

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(); 
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(); 
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(); 

    const auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    const auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    const auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, pp_prim_pair_count});

    std::vector<double> mat_S(max_prim_pair_count);
    std::vector<double> mat_T(max_prim_pair_count);
    std::vector<double> mat_V(max_prim_pair_count);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());

    // S on device

    double *d_mat_S, *d_mat_T, *d_mat_V;

    uint32_t *d_ss_first_inds_local, *d_ss_second_inds_local;
    uint32_t *d_sp_first_inds_local, *d_sp_second_inds_local;
    uint32_t *d_pp_first_inds_local, *d_pp_second_inds_local;

    hipSafe(hipMalloc(&d_mat_S, max_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_T, max_prim_pair_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_V, max_prim_pair_count * sizeof(double)));

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

    STV_matrices[0] = CDenseMatrix(naos, naos);
    STV_matrices[1] = CDenseMatrix(naos, naos);
    STV_matrices[2] = CDenseMatrix(naos, naos);

    auto& mat_overlap = STV_matrices[0];
    auto& mat_kinetic_energy = STV_matrices[1];
    auto& mat_nuclear_potential = STV_matrices[2];

    mat_overlap.zero();
    mat_kinetic_energy.zero();
    mat_nuclear_potential.zero();

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

    if (ss_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(computeOverlapAndKineticEnergySS, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local));

        hipLaunchKernelGGL(computeNuclearPotentialSS, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_ss_first_inds_local,
                                                     d_ss_second_inds_local,
                                                     static_cast<uint32_t>(ss_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
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

    if (sp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(computeOverlapAndKineticEnergySP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_sp_first_inds_local,
                                                     d_sp_second_inds_local,
                                                     static_cast<uint32_t>(sp_prim_pair_count_local));

        hipLaunchKernelGGL(computeNuclearPotentialSP, num_blocks, threads_per_block, 0, 0, d_mat_V,
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

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
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

    if (pp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x);

        hipLaunchKernelGGL(computeOverlapAndKineticEnergyPP, num_blocks, threads_per_block, 0, 0, d_mat_S,
                                                     d_mat_T,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local));

        hipLaunchKernelGGL(computeNuclearPotentialPP, num_blocks, threads_per_block, 0, 0, d_mat_V,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     d_pp_first_inds_local,
                                                     d_pp_second_inds_local,
                                                     static_cast<uint32_t>(pp_prim_pair_count_local),
                                                     d_points_info,
                                                     static_cast<uint32_t>(natoms),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(mat_S.data(), d_mat_S, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
        hipSafe(hipMemcpy(mat_T.data(), d_mat_T, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));
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

    return STV_matrices;
}

auto
computeFockOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix, CScreeningData& screening) -> CDenseMatrix
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

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

    std::string errnaos("gpu::computeCoulombFock: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

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

    timer.start("  D matrices");

    screening.updateDensityMatrix(molecule, basis, densityMatrix);

    screening.sortQD(s_prim_count, p_prim_count);

    timer.stop("  D matrices");

    timer.stop("GTO block pair prep.");

    timer.start("Coulomb prep.");

    const auto& ss_mat_Q_orig = screening.getQMatrixSS();
    const auto& sp_mat_Q_orig = screening.getQMatrixSP();
    const auto& pp_mat_Q_orig = screening.getQMatrixPP();

    const auto& ss_mat_D_orig = screening.getDMatrixSS();
    const auto& sp_mat_D_orig = screening.getDMatrixSP();
    const auto& pp_mat_D_orig = screening.getDMatrixPP();

    const auto& ss_first_inds_local = screening.get_ss_first_inds_local(); 
    const auto& sp_first_inds_local = screening.get_sp_first_inds_local(); 
    const auto& pp_first_inds_local = screening.get_pp_first_inds_local(); 

    const auto& ss_second_inds_local = screening.get_ss_second_inds_local(); 
    const auto& sp_second_inds_local = screening.get_sp_second_inds_local(); 
    const auto& pp_second_inds_local = screening.get_pp_second_inds_local(); 

    const auto& ss_mat_Q_local = screening.get_ss_mat_Q_local(); 
    const auto& sp_mat_Q_local = screening.get_sp_mat_Q_local(); 
    const auto& pp_mat_Q_local = screening.get_pp_mat_Q_local(); 

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

    const auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    const auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    const auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, pp_prim_pair_count});

    std::vector<double> mat_J(max_prim_pair_count);

    const auto ss_prim_pair_count_local = static_cast<int64_t>(ss_first_inds_local.size());
    const auto sp_prim_pair_count_local = static_cast<int64_t>(sp_first_inds_local.size());
    const auto pp_prim_pair_count_local = static_cast<int64_t>(pp_first_inds_local.size());

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
    hipSafe(hipMalloc(&d_mat_J, max_prim_pair_count * sizeof(double)));

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

    CDenseMatrix mat_Fock(naos, naos);

    mat_Fock.zero();

    timer.start("J computation");

    // compute J

    // J: (SS|SS)
    //     **

    if (ss_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SSSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SSSS");

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (SS|SP)
    //     **

    if (ss_prim_pair_count > 0 && sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SSSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SSSP");

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (SS|PP)
    //     **

    if (ss_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SSPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, ss_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SSPP");

        for (int64_t ij = 0; ij < ss_prim_pair_count_local; ij++)
        {
            const auto i = ss_first_inds_local[ij];
            const auto j = ss_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (SP|SS)
    //     **

    if (ss_prim_pair_count > 0 && sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SPSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SPSS");

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (SP|SP)
    //     **

    if (sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SPSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SPSP");

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (SP|PP)
    //     **

    if (sp_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel SPPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, sp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel SPPP");

        for (int64_t ij = 0; ij < sp_prim_pair_count_local; ij++)
        {
            const auto i = sp_first_inds_local[ij];
            const auto j = sp_second_inds_local[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (PP|SS)
    //     **

    if (ss_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel PPSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel PPSS");

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (PP|SP)
    //     **

    if (sp_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel PPSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel PPSP");

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    // J: (PP|PP)
    //     **

    if (pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("  J kernel PPPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count_local + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
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

        hipSafe(hipMemcpy(mat_J.data(), d_mat_J, pp_prim_pair_count_local * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  J kernel PPPP");

        for (int64_t ij = 0; ij < pp_prim_pair_count_local; ij++)
        {
            const auto i = pp_first_inds_local[ij];
            const auto j = pp_second_inds_local[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += mat_J[ij] * 2.0;
        }
    }

    timer.stop("J computation");

    timer.start("PreLinK Q_prime");

    // preLinK
    // J. Chem. Phys. 138, 134114 (2013)

    const auto mat_Q_full = screening.get_mat_Q_full(s_prim_count, p_prim_count);
    const auto mat_D_abs_full = screening.get_mat_D_abs_full(s_prim_count, p_prim_count);

    // TODO distribute computation of Q_prime

    double *d_matrix_A, *d_matrix_B, *d_matrix_C;

    hipSafe(hipMalloc(&d_matrix_A, mat_Q_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_B, mat_Q_full.getNumberOfElements() * sizeof(double)));
    hipSafe(hipMalloc(&d_matrix_C, mat_Q_full.getNumberOfElements() * sizeof(double)));

    hipSafe(hipMemcpy(d_matrix_A, mat_Q_full.values(), mat_Q_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_matrix_B, mat_D_abs_full.values(), mat_Q_full.getNumberOfElements() * sizeof(double), hipMemcpyHostToDevice));

    const auto all_prim_count = s_prim_count + p_prim_count * 3;

    dim3 threads_per_block(TILE_DIM, TILE_DIM);

    dim3 num_blocks((all_prim_count + threads_per_block.x - 1) / threads_per_block.x, (all_prim_count + threads_per_block.y - 1) / threads_per_block.y);

    hipLaunchKernelGGL(gpu::matmulAB, num_blocks, threads_per_block, 0, 0, d_matrix_C, d_matrix_A, d_matrix_B,
            static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count));

    hipLaunchKernelGGL(gpu::matmulAB, num_blocks, threads_per_block, 0, 0, d_matrix_B, d_matrix_C, d_matrix_A,
            static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count), static_cast<uint32_t>(all_prim_count));

    CDenseMatrix Q_prime(all_prim_count, all_prim_count);

    hipSafe(hipMemcpy(Q_prime.values(), d_matrix_B, mat_Q_full.getNumberOfElements() * sizeof(double), hipMemcpyDeviceToHost));

    // auto Q_prime_2 = denblas::multAB(mat_Q_full, denblas::multAB(mat_D_abs_full, mat_Q_full));
    // double max_diff = 0.0;
    // for (int64_t i = 0; i < mat_Q_full.getNumberOfElements(); i++)
    // {
    //     const auto diff = std::fabs(Q_prime.values()[i] - Q_prime_2.values()[i]);
    //     if (diff > max_diff) max_diff = diff;
    // }
    // printf("max_diff: %22.12e\n", max_diff);

    timer.stop("PreLinK Q_prime");

    timer.start("Exchange prep.");

    // K preparation

    const auto mat_Q_for_K_ss = screening.get_mat_Q_for_K_ss();
    const auto mat_Q_for_K_sp = screening.get_mat_Q_for_K_sp();
    const auto mat_Q_for_K_ps = screening.get_mat_Q_for_K_ps();
    const auto mat_Q_for_K_pp = screening.get_mat_Q_for_K_pp();

    const auto density_inds_for_K_ss = screening.get_density_inds_for_K_ss();
    const auto density_inds_for_K_sp = screening.get_density_inds_for_K_sp();
    const auto density_inds_for_K_ps = screening.get_density_inds_for_K_ps();
    const auto density_inds_for_K_pp = screening.get_density_inds_for_K_pp();

    const double Q_prime_thresh = 1.0e-7;

    // TODO sort by Q_prime bound

    screening.form_pair_inds_for_K(s_prim_count, p_prim_count, Q_prime, Q_prime_thresh);

    const auto pair_inds_i_for_K_ss = screening.get_pair_inds_i_for_K_ss();
    const auto pair_inds_k_for_K_ss = screening.get_pair_inds_k_for_K_ss();

    const auto pair_inds_i_for_K_sp = screening.get_pair_inds_i_for_K_sp();
    const auto pair_inds_k_for_K_sp = screening.get_pair_inds_k_for_K_sp();

    const auto pair_inds_i_for_K_pp = screening.get_pair_inds_i_for_K_pp();
    const auto pair_inds_k_for_K_pp = screening.get_pair_inds_k_for_K_pp();

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

    // K: (SS|SS)
    //     *  *

    if (ss_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_ss, 1);

        timer.start("  K kernel SSSS");

        hipLaunchKernelGGL(computeExchangeFockSSSS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SSSS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SS|SP)
    //     *  *

    if ((ss_prim_pair_count > 0) && (sp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_ss, 1);

        timer.start("  K kernel SSSP");

        hipLaunchKernelGGL(computeExchangeFockSSSP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SSSP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SP|SS)
    //     *  *

    if ((ss_prim_pair_count > 0) && (sp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_ss, 1);

        timer.start("  K kernel SPSS");

        hipLaunchKernelGGL(computeExchangeFockSPSS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SPSS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SP|SP)
    //     *  *

    if (sp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_ss, 1);

        timer.start("  K kernel SPSP");

        hipLaunchKernelGGL(computeExchangeFockSPSP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        timer.stop("  K kernel SPSP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SS|PS)
    //     *  *

    if ((ss_prim_pair_count > 0) && (sp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_sp, 1);

        timer.start("  K kernel SSPS");

        hipLaunchKernelGGL(computeExchangeFockSSPS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SSPS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);
            mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SS|PP)
    //     *  *

    if ((ss_prim_pair_count > 0) && (pp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_sp, 1);

        timer.start("  K kernel SSPP");

        hipLaunchKernelGGL(computeExchangeFockSSPP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SSPP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);
            mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SP|PS)
    //     *  *

    if (sp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_sp, 1);

        timer.start("  K kernel SPPS");

        hipLaunchKernelGGL(computeExchangeFockSPPS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_sp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel SPPS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);
            mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (SP|PP)
    //     *  *

    if ((sp_prim_pair_count > 0) && (pp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_sp, 1);

        timer.start("  K kernel SPPP");

        hipLaunchKernelGGL(computeExchangeFockSPPP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        timer.stop("  K kernel SPPP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_sp; ik++)
        {
            const auto i = pair_inds_i_for_K_sp[ik];
            const auto k = pair_inds_k_for_K_sp[ik];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);
            mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (PS|PS)
    //     *  *

    if (sp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_pp, 1);

        timer.start("  K kernel PSPS");

        hipLaunchKernelGGL(computeExchangeFockPSPS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel PSPS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (PS|PP)
    //     *  *

    if ((sp_prim_pair_count > 0) && (pp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_pp, 1);

        timer.start("  K kernel PSPP");

        hipLaunchKernelGGL(computeExchangeFockPSPP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel PSPP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (PP|PS)
    //     *  *

    if ((sp_prim_pair_count > 0) && (pp_prim_pair_count > 0))
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_pp, 1);

        timer.start("  K kernel PPPS");

        hipLaunchKernelGGL(computeExchangeFockPPPS, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_pp * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("  K kernel PPPS");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    // K: (PP|PP)
    //     *  *

    if (pp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_pp, 1);

        timer.start("  K kernel PPPP");

        hipLaunchKernelGGL(computeExchangeFockPPPP, num_blocks, threads_per_block, 0, 0, d_mat_K,
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

        timer.stop("  K kernel PPPP");

        for (int64_t ik = 0; ik < pair_inds_count_for_K_pp; ik++)
        {
            const auto i = pair_inds_i_for_K_pp[ik];
            const auto k = pair_inds_k_for_K_pp[ik];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto k_cgto = p_prim_aoinds[(k / 3) + p_prim_count * (k % 3)];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);
        }
    }

    timer.stop("K computation");

    hipSafe(hipFree(d_boys_func_table));
    hipSafe(hipFree(d_boys_func_ft));

    hipSafe(hipFree(d_s_prim_info));
    hipSafe(hipFree(d_s_prim_aoinds));

    hipSafe(hipFree(d_p_prim_info));
    hipSafe(hipFree(d_p_prim_aoinds));

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

    hipSafe(hipFree(d_matrix_A));
    hipSafe(hipFree(d_matrix_B));
    hipSafe(hipFree(d_matrix_C));

    timer.stop("Total timing");

    // std::cout << "\nTiming of ERIs on GPU\n";
    // std::cout << "-----------------------\n";
    // std::cout << timer.getSummary() << std::endl;

    return mat_Fock;
}

}  // namespace gpu
