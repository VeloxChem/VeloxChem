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

#include "GtoValuesRecD.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto get_lda_values_rec_d(const CGtoBlock&            gto_block,
                          const std::vector<double>&  grid_coords_x,
                          const std::vector<double>&  grid_coords_y,
                          const std::vector<double>&  grid_coords_z,
                          const std::vector<int>& gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);
    
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("LDA", 5 * nrows, ncols);
        
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
        
        auto submat = gto_values.sub_matrix({0, 0});
        
        // compute GTO values for D type GTOs on grid
        
        std::vector<double> buffer_xx(ncols);

        std::vector<double> buffer_xy(ncols);

        std::vector<double> buffer_xz(ncols);

        std::vector<double> buffer_yy(ncols);

        std::vector<double> buffer_yz(ncols);

        std::vector<double> buffer_zz(ncols);

        auto ptr_buffer_xx = buffer_xx.data();

        auto ptr_buffer_xy = buffer_xy.data();

        auto ptr_buffer_xz = buffer_xz.data();

        auto ptr_buffer_yy = buffer_yy.data();

        auto ptr_buffer_yz = buffer_yz.data();

        auto ptr_buffer_zz = buffer_zz.data();
        
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
                
                std::ranges::fill(buffer_xx, 0.0);
                
                std::ranges::fill(buffer_xy, 0.0);
                
                std::ranges::fill(buffer_xz, 0.0);
                
                std::ranges::fill(buffer_yy, 0.0);
                
                std::ranges::fill(buffer_yz, 0.0);
                
                std::ranges::fill(buffer_zz, 0.0);
              
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
        
                        const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        ptr_buffer_xx[k] += gr_x * gr_x * fss;

                        ptr_buffer_xy[k] += gr_x * gr_y * fss;

                        ptr_buffer_xz[k] += gr_x * gr_z * fss;

                        ptr_buffer_yy[k] += gr_y * gr_y * fss;

                        ptr_buffer_yz[k] += gr_y * gr_z * fss;

                        ptr_buffer_zz[k] += gr_z * gr_z * fss;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat, buffer_xx, -1.0, 2 * nrows + irow);

                gtoval::distribute(submat, buffer_xx, 0.5 * f2_3, 4 * nrows + irow);

                gtoval::distribute(submat, buffer_xy, f2_3, irow);

                gtoval::distribute(submat, buffer_xz, f2_3, 3 * nrows + irow);

                gtoval::distribute(submat, buffer_yy, -1.0, 2 * nrows + irow);

                gtoval::distribute(submat, buffer_yy, -0.5 * f2_3, 4 * nrows + irow);

                gtoval::distribute(submat, buffer_yz, f2_3, nrows + irow);

                gtoval::distribute(submat, buffer_zz, 2.0, 2 * nrows + irow);
    
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

auto
get_gga_values_rec_d(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", 5 * nrows, ncols);
        
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
        
        std::vector<double> buffer_0_xx(ncols);
        
        std::vector<double> buffer_0_xy(ncols);
        
        std::vector<double> buffer_0_xz(ncols);
        
        std::vector<double> buffer_0_yy(ncols);
        
        std::vector<double> buffer_0_yz(ncols);
        
        std::vector<double> buffer_0_zz(ncols);

        std::vector<double> buffer_x_xx(ncols);
        
        std::vector<double> buffer_x_xy(ncols);
        
        std::vector<double> buffer_x_xz(ncols);
        
        std::vector<double> buffer_x_yy(ncols);
        
        std::vector<double> buffer_x_yz(ncols);
        
        std::vector<double> buffer_x_zz(ncols);

        std::vector<double> buffer_y_xx(ncols);
        
        std::vector<double> buffer_y_xy(ncols);
        
        std::vector<double> buffer_y_xz(ncols);
        
        std::vector<double> buffer_y_yy(ncols);
        
        std::vector<double> buffer_y_yz(ncols);
        
        std::vector<double> buffer_y_zz(ncols);

        std::vector<double> buffer_z_xx(ncols);
        
        std::vector<double> buffer_z_xy(ncols);
        
        std::vector<double> buffer_z_xz(ncols);
        
        std::vector<double> buffer_z_yy(ncols);
        
        std::vector<double> buffer_z_yz(ncols);
        
        std::vector<double> buffer_z_zz(ncols);

        auto ptr_buffer_0_xx = buffer_0_xx.data();
        
        auto ptr_buffer_0_xy = buffer_0_xy.data();
        
        auto ptr_buffer_0_xz = buffer_0_xz.data();
        
        auto ptr_buffer_0_yy = buffer_0_yy.data();
        
        auto ptr_buffer_0_yz = buffer_0_yz.data();
        
        auto ptr_buffer_0_zz = buffer_0_zz.data();

        auto ptr_buffer_x_xx = buffer_x_xx.data();
        
        auto ptr_buffer_x_xy = buffer_x_xy.data();
        
        auto ptr_buffer_x_xz = buffer_x_xz.data();
        
        auto ptr_buffer_x_yy = buffer_x_yy.data();
        
        auto ptr_buffer_x_yz = buffer_x_yz.data();
        
        auto ptr_buffer_x_zz = buffer_x_zz.data();

        auto ptr_buffer_y_xx = buffer_y_xx.data();
        
        auto ptr_buffer_y_xy = buffer_y_xy.data();
        
        auto ptr_buffer_y_xz = buffer_y_xz.data();
        
        auto ptr_buffer_y_yy = buffer_y_yy.data();
        
        auto ptr_buffer_y_yz = buffer_y_yz.data();
        
        auto ptr_buffer_y_zz = buffer_y_zz.data();

        auto ptr_buffer_z_xx = buffer_z_xx.data();
        
        auto ptr_buffer_z_xy = buffer_z_xy.data();
        
        auto ptr_buffer_z_xz = buffer_z_xz.data();
        
        auto ptr_buffer_z_yy = buffer_z_yy.data();
        
        auto ptr_buffer_z_yz = buffer_z_yz.data();
        
        auto ptr_buffer_z_zz = buffer_z_zz.data();
        
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
                
                std::ranges::fill(buffer_0_xx, 0.0);
                
                std::ranges::fill(buffer_0_xy, 0.0);
                
                std::ranges::fill(buffer_0_xz, 0.0);
                
                std::ranges::fill(buffer_0_yy, 0.0);
                
                std::ranges::fill(buffer_0_yz, 0.0);
                
                std::ranges::fill(buffer_0_zz, 0.0);
                
                std::ranges::fill(buffer_x_xx, 0.0);
                
                std::ranges::fill(buffer_x_xy, 0.0);
                
                std::ranges::fill(buffer_x_xz, 0.0);
                
                std::ranges::fill(buffer_x_yy, 0.0);
                
                std::ranges::fill(buffer_x_yz, 0.0);
                
                std::ranges::fill(buffer_x_zz, 0.0);
                
                std::ranges::fill(buffer_y_xx, 0.0);
                
                std::ranges::fill(buffer_y_xy, 0.0);
                
                std::ranges::fill(buffer_y_xz, 0.0);
                
                std::ranges::fill(buffer_y_yy, 0.0);
                
                std::ranges::fill(buffer_y_yz, 0.0);
                
                std::ranges::fill(buffer_y_zz, 0.0);
                
                std::ranges::fill(buffer_z_xx, 0.0);
                
                std::ranges::fill(buffer_z_xy, 0.0);
                
                std::ranges::fill(buffer_z_xz, 0.0);
                
                std::ranges::fill(buffer_z_yy, 0.0);
                
                std::ranges::fill(buffer_z_yz, 0.0);
                
                std::ranges::fill(buffer_z_zz, 0.0);
                
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
                        
                        const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const auto fg0 = -2.0 * fexp;

                        const auto f0x = gr_x * f00;
                        
                        const auto f0y = gr_y * f00;
                        
                        const auto f0z = gr_z * f00;

                        ptr_buffer_0_xx[k] += f0x * gr_x;
                        
                        ptr_buffer_x_xx[k] += f0x * (2.0 + fg0 * gr_x * gr_x);
                        
                        ptr_buffer_y_xx[k] += f0x * fg0 * gr_x * gr_y;
                        
                        ptr_buffer_z_xx[k] += f0x * fg0 * gr_x * gr_z;

                        ptr_buffer_0_xy[k] += f0x * gr_y;
                        
                        ptr_buffer_x_xy[k] += f0y * (1.0 + fg0 * gr_x * gr_x);
                        
                        ptr_buffer_y_xy[k] += f0x * (1.0 + fg0 * gr_y * gr_y);
                        
                        ptr_buffer_z_xy[k] += f0x * fg0 * gr_y * gr_z;

                        ptr_buffer_0_xz[k] += f0x * gr_z;
                        
                        ptr_buffer_x_xz[k] += f0z * (1.0 + fg0 * gr_x * gr_x);
                        
                        ptr_buffer_y_xz[k] += f0x * fg0 * gr_z * gr_y;
                        
                        ptr_buffer_z_xz[k] += f0x * (1.0 + fg0 * gr_z * gr_z);

                        ptr_buffer_0_yy[k] += f0y * gr_y;
                        
                        ptr_buffer_x_yy[k] += f0y * fg0 * gr_y * gr_x;
                        
                        ptr_buffer_y_yy[k] += f0y * (2.0 + fg0 * gr_y * gr_y);
                        
                        ptr_buffer_z_yy[k] += f0y * fg0 * gr_y * gr_z;

                        ptr_buffer_0_yz[k] += f0y * gr_z;
                        
                        ptr_buffer_x_yz[k] += f0y * fg0 * gr_z * gr_x;
                        
                        ptr_buffer_y_yz[k] += f0z * (1.0 + fg0 * gr_y * gr_y);
                        
                        ptr_buffer_z_yz[k] += f0y * (1.0 + fg0 * gr_z * gr_z);

                        ptr_buffer_0_zz[k] += f0z * gr_z;
                        
                        ptr_buffer_x_zz[k] += f0z * fg0 * gr_z * gr_x;
                        
                        ptr_buffer_y_zz[k] += f0z * fg0 * gr_z * gr_y;
                        
                        ptr_buffer_z_zz[k] += f0z * (2.0 + fg0 * gr_z * gr_z);
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_xx, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_xx, 0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_xy, f2_3, irow);
                
                gtoval::distribute(submat_0, buffer_0_xz, f2_3, 3 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_yy, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_yy, -0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_yz, f2_3, nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_x, buffer_x_xx, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_xx, 0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_xy, f2_3, irow);
                
                gtoval::distribute(submat_x, buffer_x_xz, f2_3, 3 * nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_yy, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_yy, -0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_yz, f2_3, nrows + irow);
                
                gtoval::distribute(submat_x, buffer_x_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_y, buffer_y_xx, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_xx, 0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_xy, f2_3, irow);
                
                gtoval::distribute(submat_y, buffer_y_xz, f2_3, 3 * nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_yy, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_yy, -0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_yz, f2_3, nrows + irow);
                
                gtoval::distribute(submat_y, buffer_y_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_z, buffer_z_xx, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_xx, 0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_xy, f2_3, irow);
                
                gtoval::distribute(submat_z, buffer_z_xz, f2_3, 3 * nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_yy, -1.0, 2 * nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_yy, -0.5 * f2_3, 4 * nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_yz, f2_3, nrows + irow);
                
                gtoval::distribute(submat_z, buffer_z_zz, 2.0, 2 * nrows + irow);
    
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

auto
get_mgga_values_rec_d(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("MGGA", 5 * nrows, ncols);
        
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

        auto submat_xx = gto_values.sub_matrix({2, 0});
        auto submat_xy = gto_values.sub_matrix({2, 1});
        auto submat_xz = gto_values.sub_matrix({2, 2});
        auto submat_yy = gto_values.sub_matrix({2, 3});
        auto submat_yz = gto_values.sub_matrix({2, 4});
        auto submat_zz = gto_values.sub_matrix({2, 5});
        
        // compute GTO values for D type GTOs on grid
        
        std::vector<double> buffer_0_xx(ncols);
        std::vector<double> buffer_0_xy(ncols);
        std::vector<double> buffer_0_xz(ncols);
        std::vector<double> buffer_0_yy(ncols);
        std::vector<double> buffer_0_yz(ncols);
        std::vector<double> buffer_0_zz(ncols);

        std::vector<double> buffer_x_xx(ncols);
        std::vector<double> buffer_x_xy(ncols);
        std::vector<double> buffer_x_xz(ncols);
        std::vector<double> buffer_x_yy(ncols);
        std::vector<double> buffer_x_yz(ncols);
        std::vector<double> buffer_x_zz(ncols);

        std::vector<double> buffer_y_xx(ncols);
        std::vector<double> buffer_y_xy(ncols);
        std::vector<double> buffer_y_xz(ncols);
        std::vector<double> buffer_y_yy(ncols);
        std::vector<double> buffer_y_yz(ncols);
        std::vector<double> buffer_y_zz(ncols);

        std::vector<double> buffer_z_xx(ncols);
        std::vector<double> buffer_z_xy(ncols);
        std::vector<double> buffer_z_xz(ncols);
        std::vector<double> buffer_z_yy(ncols);
        std::vector<double> buffer_z_yz(ncols);
        std::vector<double> buffer_z_zz(ncols);

        std::vector<double> buffer_xx_xx(ncols);
        std::vector<double> buffer_xx_xy(ncols);
        std::vector<double> buffer_xx_xz(ncols);
        std::vector<double> buffer_xx_yy(ncols);
        std::vector<double> buffer_xx_yz(ncols);
        std::vector<double> buffer_xx_zz(ncols);

        std::vector<double> buffer_xy_xx(ncols);
        std::vector<double> buffer_xy_xy(ncols);
        std::vector<double> buffer_xy_xz(ncols);
        std::vector<double> buffer_xy_yy(ncols);
        std::vector<double> buffer_xy_yz(ncols);
        std::vector<double> buffer_xy_zz(ncols);

        std::vector<double> buffer_xz_xx(ncols);
        std::vector<double> buffer_xz_xy(ncols);
        std::vector<double> buffer_xz_xz(ncols);
        std::vector<double> buffer_xz_yy(ncols);
        std::vector<double> buffer_xz_yz(ncols);
        std::vector<double> buffer_xz_zz(ncols);

        std::vector<double> buffer_yy_xx(ncols);
        std::vector<double> buffer_yy_xy(ncols);
        std::vector<double> buffer_yy_xz(ncols);
        std::vector<double> buffer_yy_yy(ncols);
        std::vector<double> buffer_yy_yz(ncols);
        std::vector<double> buffer_yy_zz(ncols);

        std::vector<double> buffer_yz_xx(ncols);
        std::vector<double> buffer_yz_xy(ncols);
        std::vector<double> buffer_yz_xz(ncols);
        std::vector<double> buffer_yz_yy(ncols);
        std::vector<double> buffer_yz_yz(ncols);
        std::vector<double> buffer_yz_zz(ncols);

        std::vector<double> buffer_zz_xx(ncols);
        std::vector<double> buffer_zz_xy(ncols);
        std::vector<double> buffer_zz_xz(ncols);
        std::vector<double> buffer_zz_yy(ncols);
        std::vector<double> buffer_zz_yz(ncols);
        std::vector<double> buffer_zz_zz(ncols);

        auto ptr_buffer_0_xx = buffer_0_xx.data();
        auto ptr_buffer_0_xy = buffer_0_xy.data();
        auto ptr_buffer_0_xz = buffer_0_xz.data();
        auto ptr_buffer_0_yy = buffer_0_yy.data();
        auto ptr_buffer_0_yz = buffer_0_yz.data();
        auto ptr_buffer_0_zz = buffer_0_zz.data();

        auto ptr_buffer_x_xx = buffer_x_xx.data();
        auto ptr_buffer_x_xy = buffer_x_xy.data();
        auto ptr_buffer_x_xz = buffer_x_xz.data();
        auto ptr_buffer_x_yy = buffer_x_yy.data();
        auto ptr_buffer_x_yz = buffer_x_yz.data();
        auto ptr_buffer_x_zz = buffer_x_zz.data();

        auto ptr_buffer_y_xx = buffer_y_xx.data();
        auto ptr_buffer_y_xy = buffer_y_xy.data();
        auto ptr_buffer_y_xz = buffer_y_xz.data();
        auto ptr_buffer_y_yy = buffer_y_yy.data();
        auto ptr_buffer_y_yz = buffer_y_yz.data();
        auto ptr_buffer_y_zz = buffer_y_zz.data();

        auto ptr_buffer_z_xx = buffer_z_xx.data();
        auto ptr_buffer_z_xy = buffer_z_xy.data();
        auto ptr_buffer_z_xz = buffer_z_xz.data();
        auto ptr_buffer_z_yy = buffer_z_yy.data();
        auto ptr_buffer_z_yz = buffer_z_yz.data();
        auto ptr_buffer_z_zz = buffer_z_zz.data();
        
        auto ptr_buffer_xx_xx = buffer_xx_xx.data();
        auto ptr_buffer_xx_xy = buffer_xx_xy.data();
        auto ptr_buffer_xx_xz = buffer_xx_xz.data();
        auto ptr_buffer_xx_yy = buffer_xx_yy.data();
        auto ptr_buffer_xx_yz = buffer_xx_yz.data();
        auto ptr_buffer_xx_zz = buffer_xx_zz.data();

        auto ptr_buffer_xy_xx = buffer_xy_xx.data();
        auto ptr_buffer_xy_xy = buffer_xy_xy.data();
        auto ptr_buffer_xy_xz = buffer_xy_xz.data();
        auto ptr_buffer_xy_yy = buffer_xy_yy.data();
        auto ptr_buffer_xy_yz = buffer_xy_yz.data();
        auto ptr_buffer_xy_zz = buffer_xy_zz.data();

        auto ptr_buffer_xz_xx = buffer_xz_xx.data();
        auto ptr_buffer_xz_xy = buffer_xz_xy.data();
        auto ptr_buffer_xz_xz = buffer_xz_xz.data();
        auto ptr_buffer_xz_yy = buffer_xz_yy.data();
        auto ptr_buffer_xz_yz = buffer_xz_yz.data();
        auto ptr_buffer_xz_zz = buffer_xz_zz.data();
        
        auto ptr_buffer_yy_xx = buffer_yy_xx.data();
        auto ptr_buffer_yy_xy = buffer_yy_xy.data();
        auto ptr_buffer_yy_xz = buffer_yy_xz.data();
        auto ptr_buffer_yy_yy = buffer_yy_yy.data();
        auto ptr_buffer_yy_yz = buffer_yy_yz.data();
        auto ptr_buffer_yy_zz = buffer_yy_zz.data();

        auto ptr_buffer_yz_xx = buffer_yz_xx.data();
        auto ptr_buffer_yz_xy = buffer_yz_xy.data();
        auto ptr_buffer_yz_xz = buffer_yz_xz.data();
        auto ptr_buffer_yz_yy = buffer_yz_yy.data();
        auto ptr_buffer_yz_yz = buffer_yz_yz.data();
        auto ptr_buffer_yz_zz = buffer_yz_zz.data();

        auto ptr_buffer_zz_xx = buffer_zz_xx.data();
        auto ptr_buffer_zz_xy = buffer_zz_xy.data();
        auto ptr_buffer_zz_xz = buffer_zz_xz.data();
        auto ptr_buffer_zz_yy = buffer_zz_yy.data();
        auto ptr_buffer_zz_yz = buffer_zz_yz.data();
        auto ptr_buffer_zz_zz = buffer_zz_zz.data();
        
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
                
                std::ranges::fill(buffer_0_xx, 0.0);
                std::ranges::fill(buffer_0_xy, 0.0);
                std::ranges::fill(buffer_0_xz, 0.0);
                std::ranges::fill(buffer_0_yy, 0.0);
                std::ranges::fill(buffer_0_yz, 0.0);
                std::ranges::fill(buffer_0_zz, 0.0);
                
                std::ranges::fill(buffer_x_xx, 0.0);
                std::ranges::fill(buffer_x_xy, 0.0);
                std::ranges::fill(buffer_x_xz, 0.0);
                std::ranges::fill(buffer_x_yy, 0.0);
                std::ranges::fill(buffer_x_yz, 0.0);
                std::ranges::fill(buffer_x_zz, 0.0);
                
                std::ranges::fill(buffer_y_xx, 0.0);
                std::ranges::fill(buffer_y_xy, 0.0);
                std::ranges::fill(buffer_y_xz, 0.0);
                std::ranges::fill(buffer_y_yy, 0.0);
                std::ranges::fill(buffer_y_yz, 0.0);
                std::ranges::fill(buffer_y_zz, 0.0);
                
                std::ranges::fill(buffer_z_xx, 0.0);
                std::ranges::fill(buffer_z_xy, 0.0);
                std::ranges::fill(buffer_z_xz, 0.0);
                std::ranges::fill(buffer_z_yy, 0.0);
                std::ranges::fill(buffer_z_yz, 0.0);
                std::ranges::fill(buffer_z_zz, 0.0);
                
                std::ranges::fill(buffer_xx_xx, 0.0);
                std::ranges::fill(buffer_xx_xy, 0.0);
                std::ranges::fill(buffer_xx_xz, 0.0);
                std::ranges::fill(buffer_xx_yy, 0.0);
                std::ranges::fill(buffer_xx_yz, 0.0);
                std::ranges::fill(buffer_xx_zz, 0.0);
                
                std::ranges::fill(buffer_xy_xx, 0.0);
                std::ranges::fill(buffer_xy_xy, 0.0);
                std::ranges::fill(buffer_xy_xz, 0.0);
                std::ranges::fill(buffer_xy_yy, 0.0);
                std::ranges::fill(buffer_xy_yz, 0.0);
                std::ranges::fill(buffer_xy_zz, 0.0);
                
                std::ranges::fill(buffer_xz_xx, 0.0);
                std::ranges::fill(buffer_xz_xy, 0.0);
                std::ranges::fill(buffer_xz_xz, 0.0);
                std::ranges::fill(buffer_xz_yy, 0.0);
                std::ranges::fill(buffer_xz_yz, 0.0);
                std::ranges::fill(buffer_xz_zz, 0.0);
                
                std::ranges::fill(buffer_yy_xx, 0.0);
                std::ranges::fill(buffer_yy_xy, 0.0);
                std::ranges::fill(buffer_yy_xz, 0.0);
                std::ranges::fill(buffer_yy_yy, 0.0);
                std::ranges::fill(buffer_yy_yz, 0.0);
                std::ranges::fill(buffer_yy_zz, 0.0);
                
                std::ranges::fill(buffer_yz_xx, 0.0);
                std::ranges::fill(buffer_yz_xy, 0.0);
                std::ranges::fill(buffer_yz_xz, 0.0);
                std::ranges::fill(buffer_yz_yy, 0.0);
                std::ranges::fill(buffer_yz_yz, 0.0);
                std::ranges::fill(buffer_yz_zz, 0.0);
                
                std::ranges::fill(buffer_zz_xx, 0.0);
                std::ranges::fill(buffer_zz_xy, 0.0);
                std::ranges::fill(buffer_zz_xz, 0.0);
                std::ranges::fill(buffer_zz_yy, 0.0);
                std::ranges::fill(buffer_zz_yz, 0.0);
                std::ranges::fill(buffer_zz_zz, 0.0);
                
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
                        
                        const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const auto fg0 = -2.0 * fexp;

                        const auto f0x = gr_x * f00;
                        
                        const auto f0y = gr_y * f00;
                        
                        const auto f0z = gr_z * f00;

                        const auto fg1 = f00 * fg0;

                        const auto fg2 = fg1 * fg0;

                        ptr_buffer_0_xx[k] += f0x * gr_x;
                        
                        ptr_buffer_x_xx[k] += f0x * (2.0 + fg0 * gr_x * gr_x);
                        ptr_buffer_y_xx[k] += f0x * fg0 * gr_x * gr_y;
                        ptr_buffer_z_xx[k] += f0x * fg0 * gr_x * gr_z;

                        ptr_buffer_xx_xx[k] += fg2 * gr_x * gr_x * gr_x * gr_x + 5.0 * fg1 * gr_x * gr_x + 2.0 * f00;
                        ptr_buffer_xy_xx[k] += fg2 * gr_x * gr_x * gr_x * gr_y + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xz_xx[k] += fg2 * gr_x * gr_x * gr_x * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yy_xx[k] += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x;
                        ptr_buffer_yz_xx[k] += fg2 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xx[k] += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x;

                        ptr_buffer_0_xy[k] += f0x * gr_y;
                        
                        ptr_buffer_x_xy[k] += f0y * (1.0 + fg0 * gr_x * gr_x);
                        ptr_buffer_y_xy[k] += f0x * (1.0 + fg0 * gr_y * gr_y);
                        ptr_buffer_z_xy[k] += f0x * fg0 * gr_y * gr_z;

                        ptr_buffer_xx_xy[k] += fg2 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xy_xy[k] += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x + fg1 * gr_y * gr_y + f00;
                        ptr_buffer_xz_xy[k] += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
                        ptr_buffer_yy_xy[k] += fg2 * gr_x * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_yz_xy[k] += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
                        ptr_buffer_zz_xy[k] += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;

                        ptr_buffer_0_xz[k] += f0x * gr_z;
                        
                        ptr_buffer_x_xz[k] += f0z * (1.0 + fg0 * gr_x * gr_x);
                        ptr_buffer_y_xz[k] += f0x * fg0 * gr_z * gr_y;
                        ptr_buffer_z_xz[k] += f0x * (1.0 + fg0 * gr_z * gr_z);

                        ptr_buffer_xx_xz[k] += fg2 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xy_xz[k] += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
                        ptr_buffer_xz_xz[k] += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x + fg1 * gr_z * gr_z + f00;
                        ptr_buffer_yy_xz[k] += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
                        ptr_buffer_yz_xz[k] += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;
                        ptr_buffer_zz_xz[k] += fg2 * gr_x * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z;

                        ptr_buffer_0_yy[k] += f0y * gr_y;
                        
                        ptr_buffer_x_yy[k] += f0y * fg0 * gr_y * gr_x;
                        ptr_buffer_y_yy[k] += f0y * (2.0 + fg0 * gr_y * gr_y);
                        ptr_buffer_z_yy[k] += f0y * fg0 * gr_y * gr_z;

                        ptr_buffer_xx_yy[k] += fg2 * gr_x * gr_x * gr_y * gr_y + fg1 * gr_y * gr_y;
                        ptr_buffer_xy_yy[k] += fg2 * gr_x * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xz_yy[k] += fg2 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_yy_yy[k] += fg2 * gr_y * gr_y * gr_y * gr_y + 5.0 * fg1 * gr_y * gr_y + 2.0 * f00;
                        ptr_buffer_yz_yy[k] += fg2 * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_zz_yy[k] += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y;

                        ptr_buffer_0_yz[k] += f0y * gr_z;
                        
                        ptr_buffer_x_yz[k] += f0y * fg0 * gr_z * gr_x;
                        ptr_buffer_y_yz[k] += f0z * (1.0 + fg0 * gr_y * gr_y);
                        ptr_buffer_z_yz[k] += f0y * (1.0 + fg0 * gr_z * gr_z);

                        ptr_buffer_xx_yz[k] += fg2 * gr_x * gr_x * gr_y * gr_z + fg1 * gr_y * gr_z;
                        ptr_buffer_xy_yz[k] += fg2 * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_z;
                        ptr_buffer_xz_yz[k] += fg2 * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y;
                        ptr_buffer_yy_yz[k] += fg2 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_yz_yz[k] += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y + fg1 * gr_z * gr_z + f00;
                        ptr_buffer_zz_yz[k] += fg2 * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z;

                        ptr_buffer_0_zz[k] += f0z * gr_z;
                        
                        ptr_buffer_x_zz[k] += f0z * fg0 * gr_z * gr_x;
                        ptr_buffer_y_zz[k] += f0z * fg0 * gr_z * gr_y;
                        ptr_buffer_z_zz[k] += f0z * (2.0 + fg0 * gr_z * gr_z);

                        ptr_buffer_xx_zz[k] += fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_z * gr_z;
                        ptr_buffer_xy_zz[k] += fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_xz_zz[k] += fg2 * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yy_zz[k] += fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_z * gr_z;
                        ptr_buffer_yz_zz[k] += fg2 * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_zz_zz[k] += fg2 * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_z * gr_z + 2.0 * f00;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_xy, f2_3, irow);
                gtoval::distribute(submat_0, buffer_0_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_0, buffer_0_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_x, buffer_x_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_xy, f2_3, irow);
                gtoval::distribute(submat_x, buffer_x_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_x, buffer_x_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_y, buffer_y_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_xy, f2_3, irow);
                gtoval::distribute(submat_y, buffer_y_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_y, buffer_y_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_z, buffer_z_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_xy, f2_3, irow);
                gtoval::distribute(submat_z, buffer_z_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_z, buffer_z_zz, 2.0, 2 * nrows + irow);
    
                gtoval::distribute(submat_xx, buffer_xx_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_xy, f2_3, irow);
                gtoval::distribute(submat_xx, buffer_xx_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_xy, buffer_xy_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_xy, f2_3, irow);
                gtoval::distribute(submat_xy, buffer_xy_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_xz, buffer_xz_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_xy, f2_3, irow);
                gtoval::distribute(submat_xz, buffer_xz_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_zz, 2.0, 2 * nrows + irow);
    
                gtoval::distribute(submat_yy, buffer_yy_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_xy, f2_3, irow);
                gtoval::distribute(submat_yy, buffer_yy_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_yz, buffer_yz_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_xy, f2_3, irow);
                gtoval::distribute(submat_yz, buffer_yz_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_zz, 2.0, 2 * nrows + irow);

                gtoval::distribute(submat_zz, buffer_zz_xx, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_xx, 0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_xy, f2_3, irow);
                gtoval::distribute(submat_zz, buffer_zz_xz, f2_3, 3 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_yy, -1.0, 2 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_yy, -0.5 * f2_3, 4 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_yz, f2_3, nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_zz, 2.0, 2 * nrows + irow);
    
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
