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
computeCoulombFock(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> CDenseMatrix
{
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

    timer.start("Boys func. prep.");

    // Boys function (tabulated for order 0-28)

    std::vector<double> boys_func_table((28 + 1) * 121 * 7);

    for (int64_t bf_order = 0; bf_order <= 28; bf_order++)
    {
        const auto bf_table = boysfunc::getBoysFuncTable(bf_order);

        auto bf_data = boys_func_table.data() + bf_order * 121 * 7;

        for (int64_t r = 0; r < 121; r++)
        {
            for (int64_t c = 0; c < 7; c++)
            {
                bf_data[r * 7 + c] = bf_table[r][c];
            }
        }
    }

    double* d_boys_func_table;

    hipSafe(hipMalloc(&d_boys_func_table, boys_func_table.size() * sizeof(double)));

    hipSafe(hipMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

    std::vector<double> boys_func_ft{1.0,        1.0 / 3.0,  1.0 / 5.0,  1.0 / 7.0,  1.0 / 9.0,  1.0 / 11.0, 1.0 / 13.0,
                                     1.0 / 15.0, 1.0 / 17.0, 1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0, 1.0 / 27.0,
                                     1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0, 1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                                     1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0, 1.0 / 51.0, 1.0 / 53.0, 1.0 / 55.0};

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

    // S-S gto block pair and S-P gto block pair

    auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    printf("ss_prim_pair_count: %ld\n", ss_prim_pair_count);
    printf("sp_prim_pair_count: %ld\n", sp_prim_pair_count);
    printf("pp_prim_pair_count: %ld\n", pp_prim_pair_count);

    CDenseMatrix ss_mat_Q_orig(s_prim_count, s_prim_count);
    CDenseMatrix sp_mat_Q_orig(s_prim_count, p_prim_count * 3);
    CDenseMatrix pp_mat_Q_orig(p_prim_count * 3, p_prim_count * 3);

    CDenseMatrix ss_mat_D_orig(s_prim_count, s_prim_count);
    CDenseMatrix sp_mat_D_orig(s_prim_count, p_prim_count * 3);
    CDenseMatrix pp_mat_D_orig(p_prim_count * 3, p_prim_count * 3);

    std::vector<std::tuple<double, int64_t, int64_t, double>> ss_mat_Q_D_sorted;
    std::vector<std::tuple<double, int64_t, int64_t, double>> sp_mat_Q_D_sorted;
    std::vector<std::tuple<double, int64_t, int64_t, double>> pp_mat_Q_D_sorted;

    auto dens_ptr = densityMatrix.alphaDensity(0);

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto i_cgto = s_prim_aoinds[i];

        // S-S gto block pair

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

            const auto S_ij = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto F0_t = boysfunc::getBoysFunction(0.0, 0, boys_func_table.data(), boys_func_ft.data());

            const auto eri_ijij = F0_t[0] * S_ij * S_ij * std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            const auto sqrt_ijij = std::sqrt(eri_ijij);

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

            ss_mat_Q_orig.row(i)[j] = sqrt_ijij;

            if (i != j) ss_mat_Q_orig.row(j)[i] = sqrt_ijij;

            ss_mat_D_orig.row(i)[j] = D_ij;

            if (i != j) ss_mat_D_orig.row(j)[i] = D_ij;

            ss_mat_Q_D_sorted.push_back(std::make_tuple(sqrt_ijij, i, j, D_ij));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            const auto a_j = p_prim_info[j + p_prim_count * 0];
            const auto c_j = p_prim_info[j + p_prim_count * 1];
            const auto x_j = p_prim_info[j + p_prim_count * 2];
            const auto y_j = p_prim_info[j + p_prim_count * 3];
            const auto z_j = p_prim_info[j + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij_0 = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto Lambda = std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            for (int64_t s = 0; s < 3; s++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * s];

                // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

                // p-1: py
                // p_0: pz
                // p+1: px

                double S_ij_1;

                if (s == 0) S_ij_1 = -(a_i / a_ij) * (y_j - y_i) * S_ij_0;
                if (s == 1) S_ij_1 = -(a_i / a_ij) * (z_j - z_i) * S_ij_0;
                if (s == 2) S_ij_1 = -(a_i / a_ij) * (x_j - x_i) * S_ij_0;

                const auto F1_t = boysfunc::getBoysFunction(0.0, 1, boys_func_table.data(), boys_func_ft.data());

                const double eri_ijij = Lambda * (S_ij_1 * S_ij_1 * F1_t[0] + S_ij_0 * S_ij_0 * F1_t[1] / (2.0 * (a_ij + a_ij)));

                const auto sqrt_ijij = std::sqrt(eri_ijij);

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                // TODO: think about the ordering of cartesian components

                sp_mat_Q_orig.row(i)[j * 3 + s] = sqrt_ijij;

                sp_mat_D_orig.row(i)[j * 3 + s] = D_ij;

                sp_mat_Q_D_sorted.push_back(std::make_tuple(sqrt_ijij, i, j * 3 + s, D_ij));
            }
        }
    }

    // P-P gto block pair

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        const auto a_i = p_prim_info[i + p_prim_count * 0];
        const auto c_i = p_prim_info[i + p_prim_count * 1];
        const auto x_i = p_prim_info[i + p_prim_count * 2];
        const auto y_i = p_prim_info[i + p_prim_count * 3];
        const auto z_i = p_prim_info[i + p_prim_count * 4];

        for (int64_t j = i; j < p_prim_count; j++)
        {
            const auto a_j = p_prim_info[j + p_prim_count * 0];
            const auto c_j = p_prim_info[j + p_prim_count * 1];
            const auto x_j = p_prim_info[j + p_prim_count * 2];
            const auto y_j = p_prim_info[j + p_prim_count * 3];
            const auto z_j = p_prim_info[j + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_i + a_j) / (MATH_CONST_PI * (a_i + a_j + a_i + a_j)));

            const auto F2_t = boysfunc::getBoysFunction(0.0, 2, boys_func_table.data(), boys_func_ft.data());

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};

            const double rPQ[3] = {0.0, 0.0, 0.0};

            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                    const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;
                    const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

                    const auto S_ij_11 = (

                            (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) +

                            (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                            ) * S_ij_00;

                    const double eri_ijij = Lambda * (

                            S_ij_11 * S_ij_11 * F2_t[0]

                            + S_ij_11 * S_ij_00 * (a_i + a_j) / (a_i + a_j + a_i + a_j) * ( -(i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F2_t[1] )

                            + S_ij_10 * S_ij_10 / (a_i + a_j + a_i + a_j) * ( 0.5 * F2_t[1] )

                            + S_ij_10 * S_ij_01 / (a_i + a_j + a_i + a_j) * ( (j_cart == i_cart ? 1.0 : 0.0) * 0.5 * F2_t[1] )

                            + S_ij_01 * S_ij_10 / (a_i + a_j + a_i + a_j) * ( (i_cart == j_cart ? 1.0 : 0.0) * 0.5 * F2_t[1] )

                            + S_ij_01 * S_ij_01 / (a_i + a_j + a_i + a_j) * ( 0.5 * F2_t[1] )

                            + S_ij_00 * S_ij_11 * (a_i + a_j) / (a_i + a_j + a_i + a_j) * ( -(i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F2_t[1] )

                            + S_ij_00 * S_ij_00 / ((a_i + a_j + a_i + a_j) * (a_i + a_j + a_i + a_j)) * (

                                  ((i_cart == j_cart ? 1.0 : 0.0) * (i_cart == j_cart ? 1.0 : 0.0) + 1.0 +
                                   (i_cart == j_cart ? 1.0 : 0.0) * (j_cart == i_cart ? 1.0 : 0.0)) * 0.25 * F2_t[2]

                            ));

                    const auto sqrt_ijij = std::sqrt(eri_ijij);

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                    // TODO: think about the ordering of cartesian components

                    pp_mat_Q_orig.row(i * 3 + i_cart)[j * 3 + j_cart] = sqrt_ijij;

                    if (i * 3 + i_cart != j * 3 + j_cart) pp_mat_Q_orig.row(j * 3 + j_cart)[i * 3 + i_cart] = sqrt_ijij;

                    pp_mat_D_orig.row(i * 3 + i_cart)[j * 3 + j_cart] = D_ij;

                    if (i * 3 + i_cart != j * 3 + j_cart) pp_mat_D_orig.row(j * 3 + j_cart)[i * 3 + i_cart] = D_ij;

                    pp_mat_Q_D_sorted.push_back(std::make_tuple(sqrt_ijij, i * 3 + i_cart, j * 3 + j_cart, D_ij));
                }
            }
        }
    }

    std::sort(ss_mat_Q_D_sorted.begin(), ss_mat_Q_D_sorted.end());
    std::reverse(ss_mat_Q_D_sorted.begin(), ss_mat_Q_D_sorted.end());

    std::sort(sp_mat_Q_D_sorted.begin(), sp_mat_Q_D_sorted.end());
    std::reverse(sp_mat_Q_D_sorted.begin(), sp_mat_Q_D_sorted.end());

    std::sort(pp_mat_Q_D_sorted.begin(), pp_mat_Q_D_sorted.end());
    std::reverse(pp_mat_Q_D_sorted.begin(), pp_mat_Q_D_sorted.end());

    timer.stop("GTO block pair prep.");

    timer.start("Coulomb prep.");

    // sorted Q, D, and indices on host

    std::vector<double> ss_mat_Q(ss_prim_pair_count);
    std::vector<double> ss_mat_D(ss_prim_pair_count);
    std::vector<double> ss_mat_J(ss_prim_pair_count);

    std::vector<uint32_t> ss_first_inds(ss_prim_pair_count);
    std::vector<uint32_t> ss_second_inds(ss_prim_pair_count);

    double ss_max_D = 0.0;

    for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto& q_ij = ss_mat_Q_D_sorted[ij];

        auto sqrt_ijij = std::get<0>(q_ij);
        auto i         = std::get<1>(q_ij);
        auto j         = std::get<2>(q_ij);
        auto D_ij      = std::get<3>(q_ij);

        ss_mat_Q[ij] = sqrt_ijij;
        ss_mat_D[ij] = D_ij;

        ss_first_inds[ij]  = static_cast<uint32_t>(i);
        ss_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > ss_max_D) ss_max_D = std::fabs(D_ij);
    }

    std::vector<double> sp_mat_Q(sp_prim_pair_count);
    std::vector<double> sp_mat_D(sp_prim_pair_count);
    std::vector<double> sp_mat_J(sp_prim_pair_count);

    std::vector<uint32_t> sp_first_inds(sp_prim_pair_count);
    std::vector<uint32_t> sp_second_inds(sp_prim_pair_count);

    double sp_max_D = 0.0;

    for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto& q_ij = sp_mat_Q_D_sorted[ij];

        auto sqrt_ijij = std::get<0>(q_ij);
        auto i         = std::get<1>(q_ij);
        auto j         = std::get<2>(q_ij);
        auto D_ij      = std::get<3>(q_ij);

        sp_mat_Q[ij] = sqrt_ijij;
        sp_mat_D[ij] = D_ij;

        sp_first_inds[ij]  = static_cast<uint32_t>(i);
        sp_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > sp_max_D) sp_max_D = std::fabs(D_ij);
    }

    std::vector<double> pp_mat_Q(pp_prim_pair_count);
    std::vector<double> pp_mat_D(pp_prim_pair_count);
    std::vector<double> pp_mat_J(pp_prim_pair_count);

    std::vector<uint32_t> pp_first_inds(pp_prim_pair_count);
    std::vector<uint32_t> pp_second_inds(pp_prim_pair_count);

    double pp_max_D = 0.0;

    for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto& q_ij = pp_mat_Q_D_sorted[ij];

        auto sqrt_ijij = std::get<0>(q_ij);
        auto i         = std::get<1>(q_ij);
        auto j         = std::get<2>(q_ij);
        auto D_ij      = std::get<3>(q_ij);

        pp_mat_Q[ij] = sqrt_ijij;
        pp_mat_D[ij] = D_ij;

        pp_first_inds[ij]  = static_cast<uint32_t>(i);
        pp_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > pp_max_D) pp_max_D = std::fabs(D_ij);
    }

    // sorted Q, D, and indices on device

    double *  d_mat_J, *d_mat_D, *d_ss_mat_Q, *d_sp_mat_Q, *d_pp_mat_Q;
    uint32_t *d_ss_first_inds, *d_ss_second_inds;
    uint32_t *d_sp_first_inds, *d_sp_second_inds;
    uint32_t *d_pp_first_inds, *d_pp_second_inds;

    auto max_prim_pair_count = std::max(ss_prim_pair_count, std::max(sp_prim_pair_count, pp_prim_pair_count));

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

    hipSafe(hipMemcpy(d_ss_mat_Q, ss_mat_Q.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_mat_Q, sp_mat_Q.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_mat_Q, pp_mat_Q.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_ss_first_inds, ss_first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_sp_first_inds, sp_first_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pp_first_inds, pp_first_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_prim_pair_count * sizeof(uint32_t), hipMemcpyHostToDevice));

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

        timer.start("J kernel SSSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     ss_max_D,
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(ss_mat_J.data(), d_mat_J, ss_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SSSS");

        for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += ss_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += ss_mat_J[ij] * 2.0;
        }
    }

    // J: (SS|SP)
    //     **

    if (ss_prim_pair_count > 0 && sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel SSSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_sp_mat_Q,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(ss_mat_J.data(), d_mat_J, ss_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SSSP");

        for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += ss_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += ss_mat_J[ij] * 2.0;
        }
    }

    // J: (SS|PP)
    //     **

    if (ss_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel SSPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSSPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_pp_mat_Q,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(ss_mat_J.data(), d_mat_J, ss_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SSPP");

        for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];
            const auto j_cgto = s_prim_aoinds[j];

            mat_Fock.row(i_cgto)[j_cgto] += ss_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += ss_mat_J[ij] * 2.0;
        }
    }

    // J: (SP|SS)
    //     **

    if (ss_prim_pair_count > 0 && sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel SPSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     ss_max_D,
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_sp_mat_Q,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(sp_mat_J.data(), d_mat_J, sp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SPSS");

        for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += sp_mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += sp_mat_J[ij] * 2.0;
        }
    }

    // J: (SP|SP)
    //     **

    if (sp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel SPSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D,
                                                     d_sp_mat_Q,
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(sp_mat_J.data(), d_mat_J, sp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SPSP");

        for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += sp_mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += sp_mat_J[ij] * 2.0;
        }
    }

    // J: (SP|PP)
    //     **

    if (sp_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel SPPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockSPPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D,
                                                     d_sp_mat_Q,
                                                     d_pp_mat_Q,
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(sp_mat_J.data(), d_mat_J, sp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel SPPP");

        for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            const auto i_cgto = s_prim_aoinds[i];

            // TODO: think about the ordering of cartesian components
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += sp_mat_J[ij] * 2.0;
            mat_Fock.row(j_cgto)[i_cgto] += sp_mat_J[ij] * 2.0;
        }
    }

    // J: (PP|SS)
    //     **

    if (ss_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, ss_mat_D.data(), ss_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel PPSS");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPSS, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     ss_max_D,
                                                     d_mat_D,
                                                     d_ss_mat_Q,
                                                     d_pp_mat_Q,
                                                     d_ss_first_inds,
                                                     d_ss_second_inds,
                                                     static_cast<uint32_t>(ss_prim_pair_count),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(pp_mat_J.data(), d_mat_J, pp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel PPSS");

        for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += pp_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += pp_mat_J[ij] * 2.0;
        }
    }

    // J: (PP|SP)
    //     **

    if (sp_prim_pair_count > 0 && pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, sp_mat_D.data(), sp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel PPSP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPSP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_s_prim_info,
                                                     static_cast<uint32_t>(s_prim_count),
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     sp_max_D,
                                                     d_mat_D,
                                                     d_sp_mat_Q,
                                                     d_pp_mat_Q,
                                                     d_sp_first_inds,
                                                     d_sp_second_inds,
                                                     static_cast<uint32_t>(sp_prim_pair_count),
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(pp_mat_J.data(), d_mat_J, pp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel PPSP");

        for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += pp_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += pp_mat_J[ij] * 2.0;
        }
    }

    // J: (PP|PP)
    //     **

    if (pp_prim_pair_count > 0)
    {
        hipSafe(hipMemcpy(d_mat_D, pp_mat_D.data(), pp_prim_pair_count * sizeof(double), hipMemcpyHostToDevice));

        timer.start("J kernel PPPP");

        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x, 1);

        hipLaunchKernelGGL(computeCoulombFockPPPP, num_blocks, threads_per_block, 0, 0, d_mat_J,
                                                     d_p_prim_info,
                                                     static_cast<uint32_t>(p_prim_count),
                                                     pp_max_D,
                                                     d_mat_D,
                                                     d_pp_mat_Q,
                                                     d_pp_first_inds,
                                                     d_pp_second_inds,
                                                     static_cast<uint32_t>(pp_prim_pair_count),
                                                     d_boys_func_table,
                                                     d_boys_func_ft);

        hipSafe(hipMemcpy(pp_mat_J.data(), d_mat_J, pp_prim_pair_count * sizeof(double), hipMemcpyDeviceToHost));

        timer.stop("J kernel PPPP");

        for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            // TODO: think about the ordering of cartesian components
            const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            mat_Fock.row(i_cgto)[j_cgto] += pp_mat_J[ij] * 2.0;

            if (i != j) mat_Fock.row(j_cgto)[i_cgto] += pp_mat_J[ij] * 2.0;
        }
    }

    timer.stop("J computation");

    /*
    // preLinK test

    CDenseMatrix mat_Q_full(s_prim_count + p_prim_count * 3, s_prim_count + p_prim_count * 3);
    CDenseMatrix mat_D_full(s_prim_count + p_prim_count * 3, s_prim_count + p_prim_count * 3);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            mat_Q_full.row(i)[j] = ss_mat_Q_orig.row(i)[j];
            mat_D_full.row(i)[j] = ss_mat_D_orig.row(i)[j];
        }
    }

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_Q_full.row(i)[s_prim_count + j] = sp_mat_Q_orig.row(i)[j];
            mat_D_full.row(i)[s_prim_count + j] = sp_mat_D_orig.row(i)[j];

            mat_Q_full.row(s_prim_count + j)[i] = sp_mat_Q_orig.row(i)[j];
            mat_D_full.row(s_prim_count + j)[i] = sp_mat_D_orig.row(i)[j];
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_Q_full.row(s_prim_count + i)[s_prim_count + j] = pp_mat_Q_orig.row(i)[j];
            mat_D_full.row(s_prim_count + i)[s_prim_count + j] = pp_mat_D_orig.row(i)[j];
        }
    }

    auto Q_prime = denblas::multAB(mat_Q_full, denblas::multAB(mat_D_full, mat_Q_full));

    int64_t count_all = 0, count_8 = 0, count_10 = 0;

    for (int64_t i = 0; i < s_prim_count + p_prim_count * 3; i++)
    {
        for (int64_t j = 0; j < s_prim_count + p_prim_count * 3; j++)
        {
            auto val = Q_prime.row(i)[j];

            count_all++;
            if (std::fabs(val) > 1.0e-10) count_10++;
            if (std::fabs(val) > 1.0e-8) count_8++;
        }
    }

    printf("count_10 : %12.5f\n", 1.0-static_cast<double>(count_10) / static_cast<double>(count_all));
    printf("count_8  : %12.5f\n", 1.0-static_cast<double>(count_8 ) / static_cast<double>(count_all));
    */

    timer.start("Exchange prep.");

    // K preparation

    // Q_ss and D_ss for K

    std::vector<double>   mat_Q_for_K_ss(s_prim_count * s_prim_count);
    std::vector<uint32_t> density_inds_for_K_ss(s_prim_count * s_prim_count);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = ss_mat_Q_orig.row(i);

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
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            mat_Q_for_K_ss[i * s_prim_count + j]        = q_val;
            density_inds_for_K_ss[i * s_prim_count + j] = static_cast<uint32_t>(j_idx);
        }
    }

    std::vector<double>   mat_Q_for_K_sp(s_prim_count * p_prim_count * 3);
    std::vector<uint32_t> density_inds_for_K_sp(s_prim_count * p_prim_count * 3);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = sp_mat_Q_orig.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            mat_Q_for_K_sp[i * p_prim_count * 3 + j]        = q_val;
            density_inds_for_K_sp[i * p_prim_count * 3 + j] = static_cast<uint32_t>(j_idx);
        }
    }

    std::vector<uint32_t> pair_inds_i_for_K_ss;
    std::vector<uint32_t> pair_inds_k_for_K_ss;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = i; k < s_prim_count; k++)
        {
            // TODO use Q_prime (Q' = QDQ)
            pair_inds_i_for_K_ss.push_back(i);
            pair_inds_k_for_K_ss.push_back(k);
        }
    }

    std::vector<uint32_t> pair_inds_i_for_K_sp;
    std::vector<uint32_t> pair_inds_k_for_K_sp;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = 0; k < p_prim_count * 3; k++)
        {
            // TODO use Q_prime (Q' = QDQ)
            pair_inds_i_for_K_sp.push_back(i);
            pair_inds_k_for_K_sp.push_back(k);
        }
    }

    auto pair_inds_count_for_K_ss = static_cast<uint32_t>(pair_inds_i_for_K_ss.size());
    auto pair_inds_count_for_K_sp = static_cast<uint32_t>(pair_inds_i_for_K_sp.size());

    printf("pair_inds_count_for_K_ss: %d\n", pair_inds_count_for_K_ss);
    printf("pair_inds_count_for_K_sp: %d\n", pair_inds_count_for_K_sp);

    double*   d_mat_K;
    uint32_t *d_pair_inds_i_for_K_ss, *d_pair_inds_k_for_K_ss;
    uint32_t *d_pair_inds_i_for_K_sp, *d_pair_inds_k_for_K_sp;

    auto max_pair_inds_count = std::max(pair_inds_count_for_K_ss, pair_inds_count_for_K_sp);

    hipSafe(hipMalloc(&d_mat_K, max_pair_inds_count * sizeof(double)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_ss, pair_inds_count_for_K_ss * sizeof(uint32_t)));

    hipSafe(hipMalloc(&d_pair_inds_i_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_pair_inds_k_for_K_sp, pair_inds_count_for_K_sp * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_ss, pair_inds_i_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_ss, pair_inds_k_for_K_ss.data(), pair_inds_count_for_K_ss * sizeof(uint32_t), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_pair_inds_i_for_K_sp, pair_inds_i_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_pair_inds_k_for_K_sp, pair_inds_k_for_K_sp.data(), pair_inds_count_for_K_sp * sizeof(uint32_t), hipMemcpyHostToDevice));

    double* d_mat_D_full_AO;
    double *d_mat_Q_for_K_ss, *d_mat_Q_for_K_sp;
    uint32_t *d_density_inds_for_K_ss, *d_density_inds_for_K_sp;

    hipSafe(hipMalloc(&d_mat_D_full_AO, naos * naos * sizeof(double)));

    hipSafe(hipMalloc(&d_mat_Q_for_K_ss, s_prim_count * s_prim_count * sizeof(double)));
    hipSafe(hipMalloc(&d_mat_Q_for_K_sp, s_prim_count * p_prim_count * 3 * sizeof(double)));

    hipSafe(hipMalloc(&d_density_inds_for_K_ss, s_prim_count * s_prim_count * sizeof(uint32_t)));
    hipSafe(hipMalloc(&d_density_inds_for_K_sp, s_prim_count * p_prim_count * 3 * sizeof(uint32_t)));

    hipSafe(hipMemcpy(d_mat_D_full_AO, dens_ptr, naos * naos * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_mat_Q_for_K_ss, mat_Q_for_K_ss.data(), s_prim_count * s_prim_count * sizeof(double), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_mat_Q_for_K_sp, mat_Q_for_K_sp.data(), s_prim_count * p_prim_count * 3 * sizeof(double), hipMemcpyHostToDevice));

    hipSafe(hipMemcpy(d_density_inds_for_K_ss, density_inds_for_K_ss.data(), s_prim_count * s_prim_count * sizeof(uint32_t), hipMemcpyHostToDevice));
    hipSafe(hipMemcpy(d_density_inds_for_K_sp, density_inds_for_K_sp.data(), s_prim_count * p_prim_count * 3 * sizeof(uint32_t), hipMemcpyHostToDevice));

    timer.stop("Exchange prep.");

    timer.start("K computation");

    // compute K

    // K: (SS|SS)
    //     *  *

    if (ss_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM, TILE_DIM);

        dim3 num_blocks(pair_inds_count_for_K_ss, 1);

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

        std::vector<double> mat_K(pair_inds_count_for_K_ss);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

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

        std::vector<double> mat_K(pair_inds_count_for_K_ss);

        hipSafe(hipMemcpy(mat_K.data(), d_mat_K, pair_inds_count_for_K_ss * sizeof(double), hipMemcpyDeviceToHost));

        for (int64_t ik = 0; ik < pair_inds_count_for_K_ss; ik++)
        {
            const auto i = pair_inds_i_for_K_ss[ik];
            const auto k = pair_inds_k_for_K_ss[ik];

            const auto i_cgto = s_prim_aoinds[i];
            const auto k_cgto = s_prim_aoinds[k];

            mat_Fock.row(i_cgto)[k_cgto] += mat_K[ik] * (-1.0);

            if (i != k) mat_Fock.row(k_cgto)[i_cgto] += mat_K[ik] * (-1.0);

            if (std::fabs(mat_K[ik]) > 1e-10) printf(" %3d  %3d      %22.12f\n", i, k, mat_K[ik]);
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

    hipSafe(hipFree(d_mat_K));

    hipSafe(hipFree(d_pair_inds_i_for_K_ss));
    hipSafe(hipFree(d_pair_inds_k_for_K_ss));

    hipSafe(hipFree(d_mat_D_full_AO));
    hipSafe(hipFree(d_mat_Q_for_K_ss));

    hipSafe(hipFree(d_density_inds_for_K_ss));

    timer.stop("Total timing");

    std::cout << "\nTiming of ERIs on GPU\n";
    std::cout << "-----------------------\n";
    std::cout << timer.getSummary() << std::endl;

    return mat_Fock;
}

}  // namespace gpu
