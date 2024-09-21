//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "GtoValuesRecF.hpp"

#include <cmath>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
get_lda_values_rec_f(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>& gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

    // set up GTO values storage

    const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1));

    const size_t ncols = grid_coords_x.size();

    auto gto_values = matfunc::make_matrix("LDA", 7 * nrows, ncols);

    auto submat = gto_values.sub_matrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.exponents();

    const auto gto_norms = gto_block.normalization_factors();

    const auto gto_coords = gto_block.coordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.number_of_basis_functions();

    const auto npgtos = gto_block.number_of_primitives();

    // compute GTO values for S type GTOs on grid

    std::vector<double> buffer_xxx(ncols);

    std::vector<double> buffer_xxy(ncols);

    std::vector<double> buffer_xxz(ncols);

    std::vector<double> buffer_xyy(ncols);

    std::vector<double> buffer_xyz(ncols);

    std::vector<double> buffer_xzz(ncols);

    std::vector<double> buffer_yyy(ncols);

    std::vector<double> buffer_yyz(ncols);

    std::vector<double> buffer_yzz(ncols);

    std::vector<double> buffer_zzz(ncols);

    auto ptr_buffer_xxx = buffer_xxx.data();

    auto ptr_buffer_xxy = buffer_xxy.data();

    auto ptr_buffer_xxz = buffer_xxz.data();

    auto ptr_buffer_xyy = buffer_xyy.data();

    auto ptr_buffer_xyz = buffer_xyz.data();

    auto ptr_buffer_xzz = buffer_xzz.data();

    auto ptr_buffer_yyy = buffer_yyy.data();

    auto ptr_buffer_yyz = buffer_yyz.data();

    auto ptr_buffer_yzz = buffer_yzz.data();

    auto ptr_buffer_zzz = buffer_zzz.data();

    int irow = 0;

    for (int i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto xyz = gto_coords[i].coordinates();
        
            const auto r_x = xyz[0];

            const auto r_y = xyz[1];

            const auto r_z = xyz[2];

            // compute GTO values on grid

            std::ranges::fill(buffer_xxx, 0.0);

            std::ranges::fill(buffer_xxy, 0.0);

            std::ranges::fill(buffer_xxz, 0.0);

            std::ranges::fill(buffer_xyy, 0.0);

            std::ranges::fill(buffer_xyz, 0.0);

            std::ranges::fill(buffer_xzz, 0.0);

            std::ranges::fill(buffer_yyy, 0.0);

            std::ranges::fill(buffer_yyz, 0.0);

            std::ranges::fill(buffer_yzz, 0.0);

            std::ranges::fill(buffer_zzz, 0.0);

            for (int j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int k = 0; k < ncols; k++)
                {
                    const auto gr_x = g_x[k] - r_x;

                    const auto gr_y = g_y[k] - r_y;

                    const auto gr_z = g_z[k] - r_z;

                    const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    ptr_buffer_xxx[k] += gr_x * gr_x * gr_x * fss;

                    ptr_buffer_xxy[k] += gr_x * gr_x * gr_y * fss;

                    ptr_buffer_xxz[k] += gr_x * gr_x * gr_z * fss;

                    ptr_buffer_xyy[k] += gr_x * gr_y * gr_y * fss;

                    ptr_buffer_xyz[k] += gr_x * gr_y * gr_z * fss;

                    ptr_buffer_xzz[k] += gr_x * gr_z * gr_z * fss;

                    ptr_buffer_yyy[k] += gr_y * gr_y * gr_y * fss;

                    ptr_buffer_yyz[k] += gr_y * gr_y * gr_z * fss;

                    ptr_buffer_yzz[k] += gr_y * gr_z * gr_z * fss;

                    ptr_buffer_zzz[k] += gr_z * gr_z * gr_z * fss;
                }
            }

            // distribute GTO values into submatrix

            gtoval::distribute(submat, buffer_xxx, -f3_3, 4 * nrows + irow);

            gtoval::distribute(submat, buffer_xxx, f3_5, 6 * nrows + irow);

            gtoval::distribute(submat, buffer_xxy, 3.0 * f3_5, irow);

            gtoval::distribute(submat, buffer_xxy, -f3_3, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_xxz, -3.0, 3 * nrows + irow);

            gtoval::distribute(submat, buffer_xxz, 0.5 * f3_15, 5 * nrows + irow);

            gtoval::distribute(submat, buffer_xyy, -f3_3, 4 * nrows + irow);

            gtoval::distribute(submat, buffer_xyy, -3.0 * f3_5, 6 * nrows + irow);

            gtoval::distribute(submat, buffer_xyz, f3_15, nrows + irow);

            gtoval::distribute(submat, buffer_xzz, 4.0 * f3_3, 4 * nrows + irow);

            gtoval::distribute(submat, buffer_yyy, -f3_5, irow);

            gtoval::distribute(submat, buffer_yyy, -f3_3, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_yyz, -3.0, 3 * nrows + irow);

            gtoval::distribute(submat, buffer_yyz, -0.5 * f3_15, 5 * nrows + irow);

            gtoval::distribute(submat, buffer_yzz, 4.0 * f3_3, 2 * nrows + irow);

            gtoval::distribute(submat, buffer_zzz, 2.0, 3 * nrows + irow);

            irow++;
        }
    }

    return gto_values;
}

auto
get_gga_values_rec_f(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", 7 * nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x = gto_values.sub_matrix({1, 0});
        
        auto submat_y = gto_values.sub_matrix({1, 1});
        
        auto submat_z = gto_values.sub_matrix({1, 2});
        
        // compute GTO values for D type GTOs on grid
        
        std::vector<double> buffer_0_xxx(ncols);
        std::vector<double> buffer_0_xxy(ncols);
        std::vector<double> buffer_0_xxz(ncols);
        std::vector<double> buffer_0_xyy(ncols);
        std::vector<double> buffer_0_xyz(ncols);
        std::vector<double> buffer_0_xzz(ncols);
        std::vector<double> buffer_0_yyy(ncols);
        std::vector<double> buffer_0_yyz(ncols);
        std::vector<double> buffer_0_yzz(ncols);
        std::vector<double> buffer_0_zzz(ncols);

        std::vector<double> buffer_x_xxx(ncols);
        std::vector<double> buffer_x_xxy(ncols);
        std::vector<double> buffer_x_xxz(ncols);
        std::vector<double> buffer_x_xyy(ncols);
        std::vector<double> buffer_x_xyz(ncols);
        std::vector<double> buffer_x_xzz(ncols);
        std::vector<double> buffer_x_yyy(ncols);
        std::vector<double> buffer_x_yyz(ncols);
        std::vector<double> buffer_x_yzz(ncols);
        std::vector<double> buffer_x_zzz(ncols);

        std::vector<double> buffer_y_xxx(ncols);
        std::vector<double> buffer_y_xxy(ncols);
        std::vector<double> buffer_y_xxz(ncols);
        std::vector<double> buffer_y_xyy(ncols);
        std::vector<double> buffer_y_xyz(ncols);
        std::vector<double> buffer_y_xzz(ncols);
        std::vector<double> buffer_y_yyy(ncols);
        std::vector<double> buffer_y_yyz(ncols);
        std::vector<double> buffer_y_yzz(ncols);
        std::vector<double> buffer_y_zzz(ncols);

        std::vector<double> buffer_z_xxx(ncols);
        std::vector<double> buffer_z_xxy(ncols);
        std::vector<double> buffer_z_xxz(ncols);
        std::vector<double> buffer_z_xyy(ncols);
        std::vector<double> buffer_z_xyz(ncols);
        std::vector<double> buffer_z_xzz(ncols);
        std::vector<double> buffer_z_yyy(ncols);
        std::vector<double> buffer_z_yyz(ncols);
        std::vector<double> buffer_z_yzz(ncols);
        std::vector<double> buffer_z_zzz(ncols);

        auto ptr_buffer_0_xxx = buffer_0_xxx.data();
        auto ptr_buffer_0_xxy = buffer_0_xxy.data();
        auto ptr_buffer_0_xxz = buffer_0_xxz.data();
        auto ptr_buffer_0_xyy = buffer_0_xyy.data();
        auto ptr_buffer_0_xyz = buffer_0_xyz.data();
        auto ptr_buffer_0_xzz = buffer_0_xzz.data();
        auto ptr_buffer_0_yyy = buffer_0_yyy.data();
        auto ptr_buffer_0_yyz = buffer_0_yyz.data();
        auto ptr_buffer_0_yzz = buffer_0_yzz.data();
        auto ptr_buffer_0_zzz = buffer_0_zzz.data();

        auto ptr_buffer_x_xxx = buffer_x_xxx.data();
        auto ptr_buffer_x_xxy = buffer_x_xxy.data();
        auto ptr_buffer_x_xxz = buffer_x_xxz.data();
        auto ptr_buffer_x_xyy = buffer_x_xyy.data();
        auto ptr_buffer_x_xyz = buffer_x_xyz.data();
        auto ptr_buffer_x_xzz = buffer_x_xzz.data();
        auto ptr_buffer_x_yyy = buffer_x_yyy.data();
        auto ptr_buffer_x_yyz = buffer_x_yyz.data();
        auto ptr_buffer_x_yzz = buffer_x_yzz.data();
        auto ptr_buffer_x_zzz = buffer_x_zzz.data();

        auto ptr_buffer_y_xxx = buffer_y_xxx.data();
        auto ptr_buffer_y_xxy = buffer_y_xxy.data();
        auto ptr_buffer_y_xxz = buffer_y_xxz.data();
        auto ptr_buffer_y_xyy = buffer_y_xyy.data();
        auto ptr_buffer_y_xyz = buffer_y_xyz.data();
        auto ptr_buffer_y_xzz = buffer_y_xzz.data();
        auto ptr_buffer_y_yyy = buffer_y_yyy.data();
        auto ptr_buffer_y_yyz = buffer_y_yyz.data();
        auto ptr_buffer_y_yzz = buffer_y_yzz.data();
        auto ptr_buffer_y_zzz = buffer_y_zzz.data();

        auto ptr_buffer_z_xxx = buffer_z_xxx.data();
        auto ptr_buffer_z_xxy = buffer_z_xxy.data();
        auto ptr_buffer_z_xxz = buffer_z_xxz.data();
        auto ptr_buffer_z_xyy = buffer_z_xyy.data();
        auto ptr_buffer_z_xyz = buffer_z_xyz.data();
        auto ptr_buffer_z_xzz = buffer_z_xzz.data();
        auto ptr_buffer_z_yyy = buffer_z_yyy.data();
        auto ptr_buffer_z_yyz = buffer_z_yyz.data();
        auto ptr_buffer_z_yzz = buffer_z_yzz.data();
        auto ptr_buffer_z_zzz = buffer_z_zzz.data();

        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];

                std::ranges::fill(buffer_0_xxx, 0.0);
                std::ranges::fill(buffer_0_xxy, 0.0);
                std::ranges::fill(buffer_0_xxz, 0.0);
                std::ranges::fill(buffer_0_xyy, 0.0);
                std::ranges::fill(buffer_0_xyz, 0.0);
                std::ranges::fill(buffer_0_xzz, 0.0);
                std::ranges::fill(buffer_0_yyy, 0.0);
                std::ranges::fill(buffer_0_yyz, 0.0);
                std::ranges::fill(buffer_0_yzz, 0.0);
                std::ranges::fill(buffer_0_zzz, 0.0);

                std::ranges::fill(buffer_x_xxx, 0.0);
                std::ranges::fill(buffer_x_xxy, 0.0);
                std::ranges::fill(buffer_x_xxz, 0.0);
                std::ranges::fill(buffer_x_xyy, 0.0);
                std::ranges::fill(buffer_x_xyz, 0.0);
                std::ranges::fill(buffer_x_xzz, 0.0);
                std::ranges::fill(buffer_x_yyy, 0.0);
                std::ranges::fill(buffer_x_yyz, 0.0);
                std::ranges::fill(buffer_x_yzz, 0.0);
                std::ranges::fill(buffer_x_zzz, 0.0);

                std::ranges::fill(buffer_y_xxx, 0.0);
                std::ranges::fill(buffer_y_xxy, 0.0);
                std::ranges::fill(buffer_y_xxz, 0.0);
                std::ranges::fill(buffer_y_xyy, 0.0);
                std::ranges::fill(buffer_y_xyz, 0.0);
                std::ranges::fill(buffer_y_xzz, 0.0);
                std::ranges::fill(buffer_y_yyy, 0.0);
                std::ranges::fill(buffer_y_yyz, 0.0);
                std::ranges::fill(buffer_y_yzz, 0.0);
                std::ranges::fill(buffer_y_zzz, 0.0);

                std::ranges::fill(buffer_z_xxx, 0.0);
                std::ranges::fill(buffer_z_xxy, 0.0);
                std::ranges::fill(buffer_z_xxz, 0.0);
                std::ranges::fill(buffer_z_xyy, 0.0);
                std::ranges::fill(buffer_z_xyz, 0.0);
                std::ranges::fill(buffer_z_xzz, 0.0);
                std::ranges::fill(buffer_z_yyy, 0.0);
                std::ranges::fill(buffer_z_yyz, 0.0);
                std::ranges::fill(buffer_z_yzz, 0.0);
                std::ranges::fill(buffer_z_zzz, 0.0);

                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
                        
                        const auto f0_0 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const double fg_0 = -2.0 * fexp;

                        const double fg_1 = f0_0 * fg_0;

                        ptr_buffer_0_xxx[k] += f0_0 * gr_x * gr_x * gr_x;
                        ptr_buffer_x_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_x + 3.0 * f0_0 * gr_x * gr_x;
                        ptr_buffer_y_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xxy[k] += f0_0 * gr_x * gr_x * gr_y;
                        ptr_buffer_x_xxy[k] += fg_1 * gr_x * gr_x * gr_x * gr_y + 2.0 * f0_0 * gr_x * gr_y;
                        ptr_buffer_y_xxy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y + f0_0 * gr_x * gr_x;
                        ptr_buffer_z_xxy[k] += fg_1 * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xxz[k] += f0_0 * gr_x * gr_x * gr_z;
                        ptr_buffer_x_xxz[k] += fg_1 * gr_x * gr_x * gr_x * gr_z + 2.0 * f0_0 * gr_x * gr_z;
                        ptr_buffer_y_xxz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xxz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z + f0_0 * gr_x * gr_x;

                        ptr_buffer_0_xyy[k] += f0_0 * gr_x * gr_y * gr_y;
                        ptr_buffer_x_xyy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y + f0_0 * gr_y * gr_y;
                        ptr_buffer_y_xyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y + 2.0 * f0_0 * gr_x * gr_y;
                        ptr_buffer_z_xyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_0_xyz[k] += f0_0 * gr_x * gr_y * gr_z;
                        ptr_buffer_x_xyz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z + f0_0 * gr_y * gr_z;
                        ptr_buffer_y_xyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z + f0_0 * gr_x * gr_z;
                        ptr_buffer_z_xyz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z + f0_0 * gr_x * gr_y;

                        ptr_buffer_0_xzz[k] += f0_0 * gr_x * gr_z * gr_z;
                        ptr_buffer_x_xzz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z + f0_0 * gr_z * gr_z;
                        ptr_buffer_y_xzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_z_xzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z + 2.0 * f0_0 * gr_x * gr_z;

                        ptr_buffer_0_yyy[k] += f0_0 * gr_y * gr_y * gr_y;
                        ptr_buffer_x_yyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_y_yyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_y + 3.0 * f0_0 * gr_y * gr_y;
                        ptr_buffer_z_yyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yyz[k] += f0_0 * gr_y * gr_y * gr_z;
                        ptr_buffer_x_yyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_y_yyz[k] += fg_1 * gr_y * gr_y * gr_y * gr_z + 2.0 * f0_0 * gr_y * gr_z;
                        ptr_buffer_z_yyz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z + f0_0 * gr_y * gr_y;

                        ptr_buffer_0_yzz[k] += f0_0 * gr_y * gr_z * gr_z;
                        ptr_buffer_x_yzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_y_yzz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z + f0_0 * gr_z * gr_z;
                        ptr_buffer_z_yzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z + 2.0 * f0_0 * gr_y * gr_z;

                        ptr_buffer_0_zzz[k] += f0_0 * gr_z * gr_z * gr_z;
                        ptr_buffer_x_zzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_y_zzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_z_zzz[k] += fg_1 * gr_z * gr_z * gr_z * gr_z + 3.0 * f0_0 * gr_z * gr_z;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_xxx, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xxx, f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xxy, 3.0 * f3_5, irow);
                gtoval::distribute(submat_0, buffer_0_xxy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xxz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xxz, 0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xyy, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xyy, -3.0 * f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xyz, f3_15, nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xzz, 4.0 * f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yyy, -f3_5, irow);
                gtoval::distribute(submat_0, buffer_0_yyy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yyz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yyz, -0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yzz, 4.0 * f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_zzz, 2.0, 3 * nrows + irow);

                gtoval::distribute(submat_x, buffer_x_xxx, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xxx, f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xxy, 3.0 * f3_5, irow);
                gtoval::distribute(submat_x, buffer_x_xxy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xxz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xxz, 0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xyy, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xyy, -3.0 * f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xyz, f3_15, nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xzz, 4.0 * f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yyy, -f3_5, irow);
                gtoval::distribute(submat_x, buffer_x_yyy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yyz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yyz, -0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yzz, 4.0 * f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_zzz, 2.0, 3 * nrows + irow);

                gtoval::distribute(submat_y, buffer_y_xxx, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xxx, f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xxy, 3.0 * f3_5, irow);
                gtoval::distribute(submat_y, buffer_y_xxy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xxz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xxz, 0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xyy, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xyy, -3.0 * f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xyz, f3_15, nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xzz, 4.0 * f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yyy, -f3_5, irow);
                gtoval::distribute(submat_y, buffer_y_yyy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yyz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yyz, -0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yzz, 4.0 * f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_zzz, 2.0, 3 * nrows + irow);

                gtoval::distribute(submat_z, buffer_z_xxx, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xxx, f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xxy, 3.0 * f3_5, irow);
                gtoval::distribute(submat_z, buffer_z_xxy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xxz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xxz, 0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xyy, -f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xyy, -3.0 * f3_5, 6 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xyz, f3_15, nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xzz, 4.0 * f3_3, 4 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yyy, -f3_5, irow);
                gtoval::distribute(submat_z, buffer_z_yyy, -f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yyz, -3.0, 3 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yyz, -0.5 * f3_15, 5 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yzz, 4.0 * f3_3, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_zzz, 2.0, 3 * nrows + irow);


                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

}  // namespace gtoval
