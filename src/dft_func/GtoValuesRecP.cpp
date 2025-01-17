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

#include "GtoValuesRecP.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto get_lda_values_rec_p(const CGtoBlock&            gto_block,
                          const std::vector<double>&  grid_coords_x,
                          const std::vector<double>&  grid_coords_y,
                          const std::vector<double>&  grid_coords_z,
                          const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("LDA", 3 * nrows, ncols);
        
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
        
        // compute GTO values for P type GTOs on grid
        
        std::vector<double> buffer_x(ncols);

        std::vector<double> buffer_y(ncols);

        std::vector<double> buffer_z(ncols);

        auto ptr_buffer_x = buffer_x.data();

        auto ptr_buffer_y = buffer_y.data();

        auto ptr_buffer_z = buffer_z.data();
        
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
                
                std::ranges::fill(buffer_x, 0.0);
                
                std::ranges::fill(buffer_y, 0.0);
                
                std::ranges::fill(buffer_z, 0.0);
              
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

                        ptr_buffer_x[k] += gr_x * fss;

                        ptr_buffer_y[k] += gr_y * fss;

                        ptr_buffer_z[k] += gr_z * fss;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat, buffer_x, 2 * nrows + irow);

                gtoval::distribute(submat, buffer_y, irow);

                gtoval::distribute(submat, buffer_z, nrows + irow);
    
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
get_gga_values_rec_p(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", 3 * nrows, ncols);
        
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
        
        // compute GTO values for P type GTOs on grid
        
        std::vector<double> buffer_0_x(ncols);
        
        std::vector<double> buffer_0_y(ncols);
        
        std::vector<double> buffer_0_z(ncols);

        std::vector<double> buffer_x_x(ncols);
        
        std::vector<double> buffer_x_y(ncols);
        
        std::vector<double> buffer_x_z(ncols);

        std::vector<double> buffer_y_x(ncols);
        
        std::vector<double> buffer_y_y(ncols);
        
        std::vector<double> buffer_y_z(ncols);

        std::vector<double> buffer_z_x(ncols);
        
        std::vector<double> buffer_z_y(ncols);
        
        std::vector<double> buffer_z_z(ncols);

        auto ptr_buffer_0_x = buffer_0_x.data();
        
        auto ptr_buffer_0_y = buffer_0_y.data();
        
        auto ptr_buffer_0_z = buffer_0_z.data();

        auto ptr_buffer_x_x = buffer_x_x.data();
        
        auto ptr_buffer_x_y = buffer_x_y.data();
        
        auto ptr_buffer_x_z = buffer_x_z.data();

        auto ptr_buffer_y_x = buffer_y_x.data();
        
        auto ptr_buffer_y_y = buffer_y_y.data();
        
        auto ptr_buffer_y_z = buffer_y_z.data();

        auto ptr_buffer_z_x = buffer_z_x.data();
        
        auto ptr_buffer_z_y = buffer_z_y.data();
        
        auto ptr_buffer_z_z = buffer_z_z.data();
        
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
                
                std::ranges::fill(buffer_0_x, 0.0);
                
                std::ranges::fill(buffer_0_y, 0.0);
                
                std::ranges::fill(buffer_0_z, 0.0);
                
                std::ranges::fill(buffer_x_x, 0.0);
                
                std::ranges::fill(buffer_x_y, 0.0);
                
                std::ranges::fill(buffer_x_z, 0.0);
                
                std::ranges::fill(buffer_y_x, 0.0);
                
                std::ranges::fill(buffer_y_y, 0.0);
                
                std::ranges::fill(buffer_y_z, 0.0);
                
                std::ranges::fill(buffer_z_x, 0.0);
                
                std::ranges::fill(buffer_z_y, 0.0);
                
                std::ranges::fill(buffer_z_z, 0.0);
              
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

                        ptr_buffer_0_x[k] += f00 * gr_x;
                        
                        ptr_buffer_x_x[k] += f00 * (1.0 + gr_x * gr_x * fg0);
                        
                        ptr_buffer_y_x[k] += f00 * gr_x * gr_y * fg0;
                        
                        ptr_buffer_z_x[k] += f00 * gr_x * gr_z * fg0;

                        ptr_buffer_0_y[k] += f00 * gr_y;
 
                        ptr_buffer_x_y[k] += f00 * gr_y * gr_x * fg0;
                        
                        ptr_buffer_y_y[k] += f00 * (1.0 + gr_y * gr_y * fg0);
                        
                        ptr_buffer_z_y[k] += f00 * gr_y * gr_z * fg0;

                        ptr_buffer_0_z[k] += f00 * gr_z;
                        
                        ptr_buffer_x_z[k] += f00 * gr_z * gr_x * fg0;
                        
                        ptr_buffer_y_z[k] += f00 * gr_z * gr_y * fg0;
                        
                        ptr_buffer_z_z[k] += f00 * (1.0 + gr_z * gr_z * fg0);
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_x, 2 * nrows + irow);
                
                gtoval::distribute(submat_0, buffer_0_y, irow);
                
                gtoval::distribute(submat_0, buffer_0_z, nrows + irow);

                gtoval::distribute(submat_x, buffer_x_x, 2 * nrows + irow);

                gtoval::distribute(submat_x, buffer_x_y, irow);

                gtoval::distribute(submat_x, buffer_x_z, nrows + irow);

                gtoval::distribute(submat_y, buffer_y_x, 2 * nrows + irow);

                gtoval::distribute(submat_y, buffer_y_y, irow);

                gtoval::distribute(submat_y, buffer_y_z, nrows + irow);

                gtoval::distribute(submat_z, buffer_z_x, 2 * nrows + irow);

                gtoval::distribute(submat_z, buffer_z_y, irow);

                gtoval::distribute(submat_z, buffer_z_z, nrows + irow);
    
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
get_mgga_values_rec_p(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("MGGA", 3 * nrows, ncols);
        
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
        
        // compute GTO values for P type GTOs on grid
        
        std::vector<double> buffer_0_x(ncols);
        std::vector<double> buffer_0_y(ncols);
        std::vector<double> buffer_0_z(ncols);

        std::vector<double> buffer_x_x(ncols);
        std::vector<double> buffer_x_y(ncols);
        std::vector<double> buffer_x_z(ncols);

        std::vector<double> buffer_y_x(ncols);
        std::vector<double> buffer_y_y(ncols);
        std::vector<double> buffer_y_z(ncols);

        std::vector<double> buffer_z_x(ncols);
        std::vector<double> buffer_z_y(ncols);
        std::vector<double> buffer_z_z(ncols);

        std::vector<double> buffer_xx_x(ncols);
        std::vector<double> buffer_xx_y(ncols);
        std::vector<double> buffer_xx_z(ncols);

        std::vector<double> buffer_xy_x(ncols);
        std::vector<double> buffer_xy_y(ncols);
        std::vector<double> buffer_xy_z(ncols);

        std::vector<double> buffer_xz_x(ncols);
        std::vector<double> buffer_xz_y(ncols);
        std::vector<double> buffer_xz_z(ncols);

        std::vector<double> buffer_yy_x(ncols);
        std::vector<double> buffer_yy_y(ncols);
        std::vector<double> buffer_yy_z(ncols);

        std::vector<double> buffer_yz_x(ncols);
        std::vector<double> buffer_yz_y(ncols);
        std::vector<double> buffer_yz_z(ncols);

        std::vector<double> buffer_zz_x(ncols);
        std::vector<double> buffer_zz_y(ncols);
        std::vector<double> buffer_zz_z(ncols);

        auto ptr_buffer_0_x = buffer_0_x.data();
        auto ptr_buffer_0_y = buffer_0_y.data();
        auto ptr_buffer_0_z = buffer_0_z.data();

        auto ptr_buffer_x_x = buffer_x_x.data();
        auto ptr_buffer_x_y = buffer_x_y.data();
        auto ptr_buffer_x_z = buffer_x_z.data();

        auto ptr_buffer_y_x = buffer_y_x.data();
        auto ptr_buffer_y_y = buffer_y_y.data();
        auto ptr_buffer_y_z = buffer_y_z.data();

        auto ptr_buffer_z_x = buffer_z_x.data();
        auto ptr_buffer_z_y = buffer_z_y.data();
        auto ptr_buffer_z_z = buffer_z_z.data();
        
        auto ptr_buffer_xx_x = buffer_xx_x.data();
        auto ptr_buffer_xx_y = buffer_xx_y.data();
        auto ptr_buffer_xx_z = buffer_xx_z.data();

        auto ptr_buffer_xy_x = buffer_xy_x.data();
        auto ptr_buffer_xy_y = buffer_xy_y.data();
        auto ptr_buffer_xy_z = buffer_xy_z.data();

        auto ptr_buffer_xz_x = buffer_xz_x.data();
        auto ptr_buffer_xz_y = buffer_xz_y.data();
        auto ptr_buffer_xz_z = buffer_xz_z.data();
        
        auto ptr_buffer_yy_x = buffer_yy_x.data();
        auto ptr_buffer_yy_y = buffer_yy_y.data();
        auto ptr_buffer_yy_z = buffer_yy_z.data();

        auto ptr_buffer_yz_x = buffer_yz_x.data();
        auto ptr_buffer_yz_y = buffer_yz_y.data();
        auto ptr_buffer_yz_z = buffer_yz_z.data();

        auto ptr_buffer_zz_x = buffer_zz_x.data();
        auto ptr_buffer_zz_y = buffer_zz_y.data();
        auto ptr_buffer_zz_z = buffer_zz_z.data();
        
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
                
                std::ranges::fill(buffer_0_x, 0.0);
                std::ranges::fill(buffer_0_y, 0.0);
                std::ranges::fill(buffer_0_z, 0.0);
                
                std::ranges::fill(buffer_x_x, 0.0);
                std::ranges::fill(buffer_x_y, 0.0);
                std::ranges::fill(buffer_x_z, 0.0);
                
                std::ranges::fill(buffer_y_x, 0.0);
                std::ranges::fill(buffer_y_y, 0.0);
                std::ranges::fill(buffer_y_z, 0.0);
                
                std::ranges::fill(buffer_z_x, 0.0);
                std::ranges::fill(buffer_z_y, 0.0);
                std::ranges::fill(buffer_z_z, 0.0);
              
                std::ranges::fill(buffer_xx_x, 0.0);
                std::ranges::fill(buffer_xx_y, 0.0);
                std::ranges::fill(buffer_xx_z, 0.0);
                
                std::ranges::fill(buffer_xy_x, 0.0);
                std::ranges::fill(buffer_xy_y, 0.0);
                std::ranges::fill(buffer_xy_z, 0.0);
                
                std::ranges::fill(buffer_xz_x, 0.0);
                std::ranges::fill(buffer_xz_y, 0.0);
                std::ranges::fill(buffer_xz_z, 0.0);
              
                std::ranges::fill(buffer_yy_x, 0.0);
                std::ranges::fill(buffer_yy_y, 0.0);
                std::ranges::fill(buffer_yy_z, 0.0);
                
                std::ranges::fill(buffer_yz_x, 0.0);
                std::ranges::fill(buffer_yz_y, 0.0);
                std::ranges::fill(buffer_yz_z, 0.0);

                std::ranges::fill(buffer_zz_x, 0.0);
                std::ranges::fill(buffer_zz_y, 0.0);
                std::ranges::fill(buffer_zz_z, 0.0);
              
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

                        const auto fg1 = f00 * fg0;

                        const auto fg2 = fg1 * fg0;

                        ptr_buffer_0_x[k] += f00 * gr_x;
                        
                        ptr_buffer_x_x[k] += f00 * (1.0 + gr_x * gr_x * fg0);
                        ptr_buffer_y_x[k] += f00 * gr_x * gr_y * fg0;
                        ptr_buffer_z_x[k] += f00 * gr_x * gr_z * fg0;

                        ptr_buffer_xx_x[k] += fg2 * gr_x * gr_x * gr_x + 3.0 * fg1 * gr_x;
                        ptr_buffer_xy_x[k] += fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
                        ptr_buffer_xz_x[k] += fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
                        ptr_buffer_yy_x[k] += fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
                        ptr_buffer_yz_x[k] += fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_x[k] += fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;

                        ptr_buffer_0_y[k] += f00 * gr_y;
 
                        ptr_buffer_x_y[k] += f00 * gr_y * gr_x * fg0;
                        ptr_buffer_y_y[k] += f00 * (1.0 + gr_y * gr_y * fg0);
                        ptr_buffer_z_y[k] += f00 * gr_y * gr_z * fg0;

                        ptr_buffer_xx_y[k] +=  + fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
                        ptr_buffer_xy_y[k] +=  + fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
                        ptr_buffer_xz_y[k] +=  + fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_y[k] +=  + fg2 * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_y;
                        ptr_buffer_yz_y[k] +=  + fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
                        ptr_buffer_zz_y[k] +=  + fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;

                        ptr_buffer_0_z[k] += f00 * gr_z;
                        
                        ptr_buffer_x_z[k] += f00 * gr_z * gr_x * fg0;
                        ptr_buffer_y_z[k] += f00 * gr_z * gr_y * fg0;
                        ptr_buffer_z_z[k] += f00 * (1.0 + gr_z * gr_z * fg0);

                        ptr_buffer_xx_z[k] +=  + fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
                        ptr_buffer_xy_z[k] +=  + fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_z[k] +=  + fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;
                        ptr_buffer_yy_z[k] +=  + fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
                        ptr_buffer_yz_z[k] +=  + fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;
                        ptr_buffer_zz_z[k] +=  + fg2 * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_z;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_x, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_y, irow);
                gtoval::distribute(submat_0, buffer_0_z, nrows + irow);

                gtoval::distribute(submat_x, buffer_x_x, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_y, irow);
                gtoval::distribute(submat_x, buffer_x_z, nrows + irow);

                gtoval::distribute(submat_y, buffer_y_x, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_y, irow);
                gtoval::distribute(submat_y, buffer_y_z, nrows + irow);

                gtoval::distribute(submat_z, buffer_z_x, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_y, irow);
                gtoval::distribute(submat_z, buffer_z_z, nrows + irow);
    
                gtoval::distribute(submat_xx, buffer_xx_x, 2 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_y, irow);
                gtoval::distribute(submat_xx, buffer_xx_z, nrows + irow);

                gtoval::distribute(submat_xy, buffer_xy_x, 2 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_y, irow);
                gtoval::distribute(submat_xy, buffer_xy_z, nrows + irow);

                gtoval::distribute(submat_xz, buffer_xz_x, 2 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_y, irow);
                gtoval::distribute(submat_xz, buffer_xz_z, nrows + irow);
    
                gtoval::distribute(submat_yy, buffer_yy_x, 2 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_y, irow);
                gtoval::distribute(submat_yy, buffer_yy_z, nrows + irow);

                gtoval::distribute(submat_yz, buffer_yz_x, 2 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_y, irow);
                gtoval::distribute(submat_yz, buffer_yz_z, nrows + irow);

                gtoval::distribute(submat_zz, buffer_zz_x, 2 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_y, irow);
                gtoval::distribute(submat_zz, buffer_zz_z, nrows + irow);
    
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
get_3rd_order_values_rec_p(const CGtoBlock&            gto_block,
                           const std::vector<double>&  grid_coords_x,
                           const std::vector<double>&  grid_coords_y,
                           const std::vector<double>&  grid_coords_z,
                           const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("3RD_ORDER", 3 * nrows, ncols);
        
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

        auto submat_xxx = gto_values.sub_matrix({3, 0});
        auto submat_xxy = gto_values.sub_matrix({3, 1});
        auto submat_xxz = gto_values.sub_matrix({3, 2});
        auto submat_xyy = gto_values.sub_matrix({3, 3});
        auto submat_xyz = gto_values.sub_matrix({3, 4});
        auto submat_xzz = gto_values.sub_matrix({3, 5});
        auto submat_yyy = gto_values.sub_matrix({3, 6});
        auto submat_yyz = gto_values.sub_matrix({3, 7});
        auto submat_yzz = gto_values.sub_matrix({3, 8});
        auto submat_zzz = gto_values.sub_matrix({3, 9});
        
        // compute GTO values for P type GTOs on grid
        
        std::vector<double> buffer_0_x(ncols);
        std::vector<double> buffer_0_y(ncols);
        std::vector<double> buffer_0_z(ncols);

        std::vector<double> buffer_x_x(ncols);
        std::vector<double> buffer_x_y(ncols);
        std::vector<double> buffer_x_z(ncols);

        std::vector<double> buffer_y_x(ncols);
        std::vector<double> buffer_y_y(ncols);
        std::vector<double> buffer_y_z(ncols);

        std::vector<double> buffer_z_x(ncols);
        std::vector<double> buffer_z_y(ncols);
        std::vector<double> buffer_z_z(ncols);

        std::vector<double> buffer_xx_x(ncols);
        std::vector<double> buffer_xx_y(ncols);
        std::vector<double> buffer_xx_z(ncols);

        std::vector<double> buffer_xy_x(ncols);
        std::vector<double> buffer_xy_y(ncols);
        std::vector<double> buffer_xy_z(ncols);

        std::vector<double> buffer_xz_x(ncols);
        std::vector<double> buffer_xz_y(ncols);
        std::vector<double> buffer_xz_z(ncols);

        std::vector<double> buffer_yy_x(ncols);
        std::vector<double> buffer_yy_y(ncols);
        std::vector<double> buffer_yy_z(ncols);

        std::vector<double> buffer_yz_x(ncols);
        std::vector<double> buffer_yz_y(ncols);
        std::vector<double> buffer_yz_z(ncols);

        std::vector<double> buffer_zz_x(ncols);
        std::vector<double> buffer_zz_y(ncols);
        std::vector<double> buffer_zz_z(ncols);

        std::vector<double> buffer_xxx_x(ncols);
        std::vector<double> buffer_xxx_y(ncols);
        std::vector<double> buffer_xxx_z(ncols);

        std::vector<double> buffer_xxy_x(ncols);
        std::vector<double> buffer_xxy_y(ncols);
        std::vector<double> buffer_xxy_z(ncols);

        std::vector<double> buffer_xxz_x(ncols);
        std::vector<double> buffer_xxz_y(ncols);
        std::vector<double> buffer_xxz_z(ncols);

        std::vector<double> buffer_xyy_x(ncols);
        std::vector<double> buffer_xyy_y(ncols);
        std::vector<double> buffer_xyy_z(ncols);

        std::vector<double> buffer_xyz_x(ncols);
        std::vector<double> buffer_xyz_y(ncols);
        std::vector<double> buffer_xyz_z(ncols);

        std::vector<double> buffer_xzz_x(ncols);
        std::vector<double> buffer_xzz_y(ncols);
        std::vector<double> buffer_xzz_z(ncols);

        std::vector<double> buffer_yyy_x(ncols);
        std::vector<double> buffer_yyy_y(ncols);
        std::vector<double> buffer_yyy_z(ncols);

        std::vector<double> buffer_yyz_x(ncols);
        std::vector<double> buffer_yyz_y(ncols);
        std::vector<double> buffer_yyz_z(ncols);

        std::vector<double> buffer_yzz_x(ncols);
        std::vector<double> buffer_yzz_y(ncols);
        std::vector<double> buffer_yzz_z(ncols);

        std::vector<double> buffer_zzz_x(ncols);
        std::vector<double> buffer_zzz_y(ncols);
        std::vector<double> buffer_zzz_z(ncols);

        auto ptr_buffer_0_x = buffer_0_x.data();
        auto ptr_buffer_0_y = buffer_0_y.data();
        auto ptr_buffer_0_z = buffer_0_z.data();

        auto ptr_buffer_x_x = buffer_x_x.data();
        auto ptr_buffer_x_y = buffer_x_y.data();
        auto ptr_buffer_x_z = buffer_x_z.data();

        auto ptr_buffer_y_x = buffer_y_x.data();
        auto ptr_buffer_y_y = buffer_y_y.data();
        auto ptr_buffer_y_z = buffer_y_z.data();

        auto ptr_buffer_z_x = buffer_z_x.data();
        auto ptr_buffer_z_y = buffer_z_y.data();
        auto ptr_buffer_z_z = buffer_z_z.data();
        
        auto ptr_buffer_xx_x = buffer_xx_x.data();
        auto ptr_buffer_xx_y = buffer_xx_y.data();
        auto ptr_buffer_xx_z = buffer_xx_z.data();

        auto ptr_buffer_xy_x = buffer_xy_x.data();
        auto ptr_buffer_xy_y = buffer_xy_y.data();
        auto ptr_buffer_xy_z = buffer_xy_z.data();

        auto ptr_buffer_xz_x = buffer_xz_x.data();
        auto ptr_buffer_xz_y = buffer_xz_y.data();
        auto ptr_buffer_xz_z = buffer_xz_z.data();
        
        auto ptr_buffer_yy_x = buffer_yy_x.data();
        auto ptr_buffer_yy_y = buffer_yy_y.data();
        auto ptr_buffer_yy_z = buffer_yy_z.data();

        auto ptr_buffer_yz_x = buffer_yz_x.data();
        auto ptr_buffer_yz_y = buffer_yz_y.data();
        auto ptr_buffer_yz_z = buffer_yz_z.data();

        auto ptr_buffer_zz_x = buffer_zz_x.data();
        auto ptr_buffer_zz_y = buffer_zz_y.data();
        auto ptr_buffer_zz_z = buffer_zz_z.data();
        
        auto ptr_buffer_xxx_x = buffer_xxx_x.data();
        auto ptr_buffer_xxx_y = buffer_xxx_y.data();
        auto ptr_buffer_xxx_z = buffer_xxx_z.data();

        auto ptr_buffer_xxy_x = buffer_xxy_x.data();
        auto ptr_buffer_xxy_y = buffer_xxy_y.data();
        auto ptr_buffer_xxy_z = buffer_xxy_z.data();

        auto ptr_buffer_xxz_x = buffer_xxz_x.data();
        auto ptr_buffer_xxz_y = buffer_xxz_y.data();
        auto ptr_buffer_xxz_z = buffer_xxz_z.data();
        
        auto ptr_buffer_xyy_x = buffer_xyy_x.data();
        auto ptr_buffer_xyy_y = buffer_xyy_y.data();
        auto ptr_buffer_xyy_z = buffer_xyy_z.data();

        auto ptr_buffer_xyz_x = buffer_xyz_x.data();
        auto ptr_buffer_xyz_y = buffer_xyz_y.data();
        auto ptr_buffer_xyz_z = buffer_xyz_z.data();

        auto ptr_buffer_xzz_x = buffer_xzz_x.data();
        auto ptr_buffer_xzz_y = buffer_xzz_y.data();
        auto ptr_buffer_xzz_z = buffer_xzz_z.data();
        
        auto ptr_buffer_yyy_x = buffer_yyy_x.data();
        auto ptr_buffer_yyy_y = buffer_yyy_y.data();
        auto ptr_buffer_yyy_z = buffer_yyy_z.data();

        auto ptr_buffer_yyz_x = buffer_yyz_x.data();
        auto ptr_buffer_yyz_y = buffer_yyz_y.data();
        auto ptr_buffer_yyz_z = buffer_yyz_z.data();

        auto ptr_buffer_yzz_x = buffer_yzz_x.data();
        auto ptr_buffer_yzz_y = buffer_yzz_y.data();
        auto ptr_buffer_yzz_z = buffer_yzz_z.data();

        auto ptr_buffer_zzz_x = buffer_zzz_x.data();
        auto ptr_buffer_zzz_y = buffer_zzz_y.data();
        auto ptr_buffer_zzz_z = buffer_zzz_z.data();
        
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
                
                std::ranges::fill(buffer_0_x, 0.0);
                std::ranges::fill(buffer_0_y, 0.0);
                std::ranges::fill(buffer_0_z, 0.0);
                
                std::ranges::fill(buffer_x_x, 0.0);
                std::ranges::fill(buffer_x_y, 0.0);
                std::ranges::fill(buffer_x_z, 0.0);
                
                std::ranges::fill(buffer_y_x, 0.0);
                std::ranges::fill(buffer_y_y, 0.0);
                std::ranges::fill(buffer_y_z, 0.0);
                
                std::ranges::fill(buffer_z_x, 0.0);
                std::ranges::fill(buffer_z_y, 0.0);
                std::ranges::fill(buffer_z_z, 0.0);
              
                std::ranges::fill(buffer_xx_x, 0.0);
                std::ranges::fill(buffer_xx_y, 0.0);
                std::ranges::fill(buffer_xx_z, 0.0);
                
                std::ranges::fill(buffer_xy_x, 0.0);
                std::ranges::fill(buffer_xy_y, 0.0);
                std::ranges::fill(buffer_xy_z, 0.0);
                
                std::ranges::fill(buffer_xz_x, 0.0);
                std::ranges::fill(buffer_xz_y, 0.0);
                std::ranges::fill(buffer_xz_z, 0.0);
              
                std::ranges::fill(buffer_yy_x, 0.0);
                std::ranges::fill(buffer_yy_y, 0.0);
                std::ranges::fill(buffer_yy_z, 0.0);
                
                std::ranges::fill(buffer_yz_x, 0.0);
                std::ranges::fill(buffer_yz_y, 0.0);
                std::ranges::fill(buffer_yz_z, 0.0);

                std::ranges::fill(buffer_zz_x, 0.0);
                std::ranges::fill(buffer_zz_y, 0.0);
                std::ranges::fill(buffer_zz_z, 0.0);
              
                std::ranges::fill(buffer_xxx_x, 0.0);
                std::ranges::fill(buffer_xxx_y, 0.0);
                std::ranges::fill(buffer_xxx_z, 0.0);
                
                std::ranges::fill(buffer_xxy_x, 0.0);
                std::ranges::fill(buffer_xxy_y, 0.0);
                std::ranges::fill(buffer_xxy_z, 0.0);
                
                std::ranges::fill(buffer_xxz_x, 0.0);
                std::ranges::fill(buffer_xxz_y, 0.0);
                std::ranges::fill(buffer_xxz_z, 0.0);
              
                std::ranges::fill(buffer_xyy_x, 0.0);
                std::ranges::fill(buffer_xyy_y, 0.0);
                std::ranges::fill(buffer_xyy_z, 0.0);
                
                std::ranges::fill(buffer_xyz_x, 0.0);
                std::ranges::fill(buffer_xyz_y, 0.0);
                std::ranges::fill(buffer_xyz_z, 0.0);

                std::ranges::fill(buffer_xzz_x, 0.0);
                std::ranges::fill(buffer_xzz_y, 0.0);
                std::ranges::fill(buffer_xzz_z, 0.0);
              
                std::ranges::fill(buffer_yyy_x, 0.0);
                std::ranges::fill(buffer_yyy_y, 0.0);
                std::ranges::fill(buffer_yyy_z, 0.0);
                
                std::ranges::fill(buffer_yyz_x, 0.0);
                std::ranges::fill(buffer_yyz_y, 0.0);
                std::ranges::fill(buffer_yyz_z, 0.0);

                std::ranges::fill(buffer_yzz_x, 0.0);
                std::ranges::fill(buffer_yzz_y, 0.0);
                std::ranges::fill(buffer_yzz_z, 0.0);

                std::ranges::fill(buffer_zzz_x, 0.0);
                std::ranges::fill(buffer_zzz_y, 0.0);
                std::ranges::fill(buffer_zzz_z, 0.0);
              
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

                        const auto fg1 = f00 * fg0;

                        const auto fg2 = fg1 * fg0;

                        const auto fg3 = fg2 * fg0;

                        ptr_buffer_0_x[k] += f00 * gr_x;
                        
                        ptr_buffer_x_x[k] += f00 * (1.0 + gr_x * gr_x * fg0);
                        ptr_buffer_y_x[k] += f00 * gr_x * gr_y * fg0;
                        ptr_buffer_z_x[k] += f00 * gr_x * gr_z * fg0;

                        ptr_buffer_xx_x[k] += fg2 * gr_x * gr_x * gr_x + 3.0 * fg1 * gr_x;
                        ptr_buffer_xy_x[k] += fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
                        ptr_buffer_xz_x[k] += fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
                        ptr_buffer_yy_x[k] += fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
                        ptr_buffer_yz_x[k] += fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_x[k] += fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;

                        ptr_buffer_xxx_x[k] += fg3 * gr_x * gr_x * gr_x * gr_x + 6.0 * fg2 * gr_x * gr_x + 3.0 * fg1;
                        ptr_buffer_xxy_x[k] += fg3 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg2 * gr_x * gr_y;
                        ptr_buffer_xxz_x[k] += fg3 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg2 * gr_x * gr_z;
                        ptr_buffer_xyy_x[k] += fg3 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x + fg2 * gr_y * gr_y + fg1;
                        ptr_buffer_xyz_x[k] += fg3 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z;
                        ptr_buffer_xzz_x[k] += fg3 * gr_x * gr_x * gr_z * gr_z + fg2 * gr_x * gr_x + fg2 * gr_z * gr_z + fg1;
                        ptr_buffer_yyy_x[k] += fg3 * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_y;
                        ptr_buffer_yyz_x[k] += fg3 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z;
                        ptr_buffer_yzz_x[k] += fg3 * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y;
                        ptr_buffer_zzz_x[k] += fg3 * gr_x * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_z;

                        ptr_buffer_0_y[k] += f00 * gr_y;
 
                        ptr_buffer_x_y[k] += f00 * gr_y * gr_x * fg0;
                        ptr_buffer_y_y[k] += f00 * (1.0 + gr_y * gr_y * fg0);
                        ptr_buffer_z_y[k] += f00 * gr_y * gr_z * fg0;

                        ptr_buffer_xx_y[k] += fg2 * gr_x * gr_x * gr_y + fg1 * gr_y;
                        ptr_buffer_xy_y[k] += fg2 * gr_x * gr_y * gr_y + fg1 * gr_x;
                        ptr_buffer_xz_y[k] += fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_y[k] += fg2 * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_y;
                        ptr_buffer_yz_y[k] += fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
                        ptr_buffer_zz_y[k] += fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;

                        ptr_buffer_xxx_y[k] += fg3 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg2 * gr_x * gr_y;
                        ptr_buffer_xxy_y[k] += fg3 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x + fg2 * gr_y * gr_y + fg1;
                        ptr_buffer_xxz_y[k] += fg3 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z;
                        ptr_buffer_xyy_y[k] += fg3 * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_y;
                        ptr_buffer_xyz_y[k] += fg3 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z;
                        ptr_buffer_xzz_y[k] += fg3 * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y;
                        ptr_buffer_yyy_y[k] += fg3 * gr_y * gr_y * gr_y * gr_y + 6.0 * fg2 * gr_y * gr_y + 3.0 * fg1;
                        ptr_buffer_yyz_y[k] += fg3 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_y * gr_z;
                        ptr_buffer_yzz_y[k] += fg3 * gr_y * gr_y * gr_z * gr_z + fg2 * gr_y * gr_y + fg2 * gr_z * gr_z + fg1;
                        ptr_buffer_zzz_y[k] += fg3 * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_z;

                        ptr_buffer_0_z[k] += f00 * gr_z;
                        
                        ptr_buffer_x_z[k] += f00 * gr_z * gr_x * fg0;
                        ptr_buffer_y_z[k] += f00 * gr_z * gr_y * fg0;
                        ptr_buffer_z_z[k] += f00 * (1.0 + gr_z * gr_z * fg0);

                        ptr_buffer_xx_z[k] +=  + fg2 * gr_x * gr_x * gr_z + fg1 * gr_z;
                        ptr_buffer_xy_z[k] +=  + fg2 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_z[k] +=  + fg2 * gr_x * gr_z * gr_z + fg1 * gr_x;
                        ptr_buffer_yy_z[k] +=  + fg2 * gr_y * gr_y * gr_z + fg1 * gr_z;
                        ptr_buffer_yz_z[k] +=  + fg2 * gr_y * gr_z * gr_z + fg1 * gr_y;
                        ptr_buffer_zz_z[k] +=  + fg2 * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_z;

                        ptr_buffer_xxx_z[k] += fg3 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg2 * gr_x * gr_z;
                        ptr_buffer_xxy_z[k] += fg3 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z;
                        ptr_buffer_xxz_z[k] += fg3 * gr_x * gr_x * gr_z * gr_z + fg2 * gr_x * gr_x + fg2 * gr_z * gr_z + fg1;
                        ptr_buffer_xyy_z[k] += fg3 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z;
                        ptr_buffer_xyz_z[k] += fg3 * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y;
                        ptr_buffer_xzz_z[k] += fg3 * gr_x * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_z;
                        ptr_buffer_yyy_z[k] += fg3 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_y * gr_z;
                        ptr_buffer_yyz_z[k] += fg3 * gr_y * gr_y * gr_z * gr_z + fg2 * gr_y * gr_y + fg2 * gr_z * gr_z + fg1;
                        ptr_buffer_yzz_z[k] += fg3 * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_z;
                        ptr_buffer_zzz_z[k] += fg3 * gr_z * gr_z * gr_z * gr_z + 6.0 * fg2 * gr_z * gr_z + 3.0 * fg1;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0, buffer_0_x, 2 * nrows + irow);
                gtoval::distribute(submat_0, buffer_0_y, irow);
                gtoval::distribute(submat_0, buffer_0_z, nrows + irow);

                gtoval::distribute(submat_x, buffer_x_x, 2 * nrows + irow);
                gtoval::distribute(submat_x, buffer_x_y, irow);
                gtoval::distribute(submat_x, buffer_x_z, nrows + irow);

                gtoval::distribute(submat_y, buffer_y_x, 2 * nrows + irow);
                gtoval::distribute(submat_y, buffer_y_y, irow);
                gtoval::distribute(submat_y, buffer_y_z, nrows + irow);

                gtoval::distribute(submat_z, buffer_z_x, 2 * nrows + irow);
                gtoval::distribute(submat_z, buffer_z_y, irow);
                gtoval::distribute(submat_z, buffer_z_z, nrows + irow);
    
                gtoval::distribute(submat_xx, buffer_xx_x, 2 * nrows + irow);
                gtoval::distribute(submat_xx, buffer_xx_y, irow);
                gtoval::distribute(submat_xx, buffer_xx_z, nrows + irow);

                gtoval::distribute(submat_xy, buffer_xy_x, 2 * nrows + irow);
                gtoval::distribute(submat_xy, buffer_xy_y, irow);
                gtoval::distribute(submat_xy, buffer_xy_z, nrows + irow);

                gtoval::distribute(submat_xz, buffer_xz_x, 2 * nrows + irow);
                gtoval::distribute(submat_xz, buffer_xz_y, irow);
                gtoval::distribute(submat_xz, buffer_xz_z, nrows + irow);
    
                gtoval::distribute(submat_yy, buffer_yy_x, 2 * nrows + irow);
                gtoval::distribute(submat_yy, buffer_yy_y, irow);
                gtoval::distribute(submat_yy, buffer_yy_z, nrows + irow);

                gtoval::distribute(submat_yz, buffer_yz_x, 2 * nrows + irow);
                gtoval::distribute(submat_yz, buffer_yz_y, irow);
                gtoval::distribute(submat_yz, buffer_yz_z, nrows + irow);

                gtoval::distribute(submat_zz, buffer_zz_x, 2 * nrows + irow);
                gtoval::distribute(submat_zz, buffer_zz_y, irow);
                gtoval::distribute(submat_zz, buffer_zz_z, nrows + irow);
    
                gtoval::distribute(submat_xxx, buffer_xxx_x, 2 * nrows + irow);
                gtoval::distribute(submat_xxx, buffer_xxx_y, irow);
                gtoval::distribute(submat_xxx, buffer_xxx_z, nrows + irow);

                gtoval::distribute(submat_xxy, buffer_xxy_x, 2 * nrows + irow);
                gtoval::distribute(submat_xxy, buffer_xxy_y, irow);
                gtoval::distribute(submat_xxy, buffer_xxy_z, nrows + irow);

                gtoval::distribute(submat_xxz, buffer_xxz_x, 2 * nrows + irow);
                gtoval::distribute(submat_xxz, buffer_xxz_y, irow);
                gtoval::distribute(submat_xxz, buffer_xxz_z, nrows + irow);
    
                gtoval::distribute(submat_xyy, buffer_xyy_x, 2 * nrows + irow);
                gtoval::distribute(submat_xyy, buffer_xyy_y, irow);
                gtoval::distribute(submat_xyy, buffer_xyy_z, nrows + irow);

                gtoval::distribute(submat_xyz, buffer_xyz_x, 2 * nrows + irow);
                gtoval::distribute(submat_xyz, buffer_xyz_y, irow);
                gtoval::distribute(submat_xyz, buffer_xyz_z, nrows + irow);

                gtoval::distribute(submat_xzz, buffer_xzz_x, 2 * nrows + irow);
                gtoval::distribute(submat_xzz, buffer_xzz_y, irow);
                gtoval::distribute(submat_xzz, buffer_xzz_z, nrows + irow);
    
                gtoval::distribute(submat_yyy, buffer_yyy_x, 2 * nrows + irow);
                gtoval::distribute(submat_yyy, buffer_yyy_y, irow);
                gtoval::distribute(submat_yyy, buffer_yyy_z, nrows + irow);

                gtoval::distribute(submat_yyz, buffer_yyz_x, 2 * nrows + irow);
                gtoval::distribute(submat_yyz, buffer_yyz_y, irow);
                gtoval::distribute(submat_yyz, buffer_yyz_z, nrows + irow);

                gtoval::distribute(submat_yzz, buffer_yzz_x, 2 * nrows + irow);
                gtoval::distribute(submat_yzz, buffer_yzz_y, irow);
                gtoval::distribute(submat_yzz, buffer_yzz_z, nrows + irow);

                gtoval::distribute(submat_zzz, buffer_zzz_x, 2 * nrows + irow);
                gtoval::distribute(submat_zzz, buffer_zzz_y, irow);
                gtoval::distribute(submat_zzz, buffer_zzz_z, nrows + irow);
    
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
