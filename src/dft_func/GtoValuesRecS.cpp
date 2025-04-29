//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "GtoValuesRecS.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
get_lda_values_rec_s(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("LDA", nrows, ncols);
        
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
        
        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer(ncols);
        
        auto ptr_buffer = buffer.data();
        
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
                
                std::ranges::fill(buffer, 0.0);
              
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
        
                        ptr_buffer[k] += fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat, buffer, irow);
    
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
get_gga_values_rec_s(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", nrows, ncols);
        
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
        
        auto submat_0_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x_0 = gto_values.sub_matrix({1, 0});
        
        auto submat_y_0 = gto_values.sub_matrix({1, 1});
        
        auto submat_z_0 = gto_values.sub_matrix({1, 2});
        
        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer_0_0(ncols);
        
        std::vector<double> buffer_x_0(ncols);
        
        std::vector<double> buffer_y_0(ncols);
        
        std::vector<double> buffer_z_0(ncols);
        
        auto ptr_buffer_0_0 = buffer_0_0.data();
        
        auto ptr_buffer_x_0 = buffer_x_0.data();
        
        auto ptr_buffer_y_0 = buffer_y_0.data();
        
        auto ptr_buffer_z_0 = buffer_z_0.data();
        
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
                
                std::ranges::fill(buffer_0_0, 0.0);
                
                std::ranges::fill(buffer_x_0, 0.0);
                
                std::ranges::fill(buffer_y_0, 0.0);
                
                std::ranges::fill(buffer_z_0, 0.0);
              
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
                        
                        const auto g0 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
                        
                        const auto g1 = -2.0 * fexp * g0;
                        
                        ptr_buffer_0_0[k] += g0;
                        
                        ptr_buffer_x_0[k] += gr_x * g1;
                        
                        ptr_buffer_y_0[k] += gr_y * g1;
                                        
                        ptr_buffer_z_0[k] += gr_z * g1;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0_0, buffer_0_0, irow);
                
                gtoval::distribute(submat_x_0, buffer_x_0, irow);
                
                gtoval::distribute(submat_y_0, buffer_y_0, irow);
                
                gtoval::distribute(submat_z_0, buffer_z_0, irow);
    
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
get_mgga_values_rec_s(const CGtoBlock&            gto_block,
                      const std::vector<double>&  grid_coords_x,
                      const std::vector<double>&  grid_coords_y,
                      const std::vector<double>&  grid_coords_z,
                      const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("MGGA", nrows, ncols);
        
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
        
        auto submat_0_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x_0 = gto_values.sub_matrix({1, 0});
        auto submat_y_0 = gto_values.sub_matrix({1, 1});
        auto submat_z_0 = gto_values.sub_matrix({1, 2});

        auto submat_xx_0 = gto_values.sub_matrix({2, 0});
        auto submat_xy_0 = gto_values.sub_matrix({2, 1});
        auto submat_xz_0 = gto_values.sub_matrix({2, 2});
        auto submat_yy_0 = gto_values.sub_matrix({2, 3});
        auto submat_yz_0 = gto_values.sub_matrix({2, 4});
        auto submat_zz_0 = gto_values.sub_matrix({2, 5});

        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer_0_0(ncols);
        
        std::vector<double> buffer_x_0(ncols);
        std::vector<double> buffer_y_0(ncols);
        std::vector<double> buffer_z_0(ncols);

        std::vector<double> buffer_xx_0(ncols);
        std::vector<double> buffer_xy_0(ncols);
        std::vector<double> buffer_xz_0(ncols);
        std::vector<double> buffer_yy_0(ncols);
        std::vector<double> buffer_yz_0(ncols);
        std::vector<double> buffer_zz_0(ncols);
        
        auto ptr_buffer_0_0 = buffer_0_0.data();
        
        auto ptr_buffer_x_0 = buffer_x_0.data();
        auto ptr_buffer_y_0 = buffer_y_0.data();
        auto ptr_buffer_z_0 = buffer_z_0.data();

        auto ptr_buffer_xx_0 = buffer_xx_0.data();
        auto ptr_buffer_xy_0 = buffer_xy_0.data();
        auto ptr_buffer_xz_0 = buffer_xz_0.data();
        auto ptr_buffer_yy_0 = buffer_yy_0.data();
        auto ptr_buffer_yz_0 = buffer_yz_0.data();
        auto ptr_buffer_zz_0 = buffer_zz_0.data();
        
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
                
                std::ranges::fill(buffer_0_0, 0.0);
                
                std::ranges::fill(buffer_x_0, 0.0);
                std::ranges::fill(buffer_y_0, 0.0);
                std::ranges::fill(buffer_z_0, 0.0);

                std::ranges::fill(buffer_xx_0, 0.0);
                std::ranges::fill(buffer_xy_0, 0.0);
                std::ranges::fill(buffer_xz_0, 0.0);
                std::ranges::fill(buffer_yy_0, 0.0);
                std::ranges::fill(buffer_yz_0, 0.0);
                std::ranges::fill(buffer_zz_0, 0.0);
              
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

                        const auto fg_0 = -2.0 * fexp;

                        const auto fg_1 = f0_0 * fg_0;

                        const auto fg_2 = fg_1 * fg_0;
                        
                        ptr_buffer_0_0[k] += f0_0;
                        
                        ptr_buffer_x_0[k] += gr_x * fg_1;
                        ptr_buffer_y_0[k] += gr_y * fg_1;
                        ptr_buffer_z_0[k] += gr_z * fg_1;

                        ptr_buffer_xx_0[k] += fg_2 * gr_x * gr_x + fg_1;
                        ptr_buffer_xy_0[k] += fg_2 * gr_x * gr_y;
                        ptr_buffer_xz_0[k] += fg_2 * gr_x * gr_z;
                        ptr_buffer_yy_0[k] += fg_2 * gr_y * gr_y + fg_1;
                        ptr_buffer_yz_0[k] += fg_2 * gr_y * gr_z;
                        ptr_buffer_zz_0[k] += fg_2 * gr_z * gr_z + fg_1;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0_0, buffer_0_0, irow);
                
                gtoval::distribute(submat_x_0, buffer_x_0, irow);
                gtoval::distribute(submat_y_0, buffer_y_0, irow);
                gtoval::distribute(submat_z_0, buffer_z_0, irow);

                gtoval::distribute(submat_xx_0, buffer_xx_0, irow);
                gtoval::distribute(submat_xy_0, buffer_xy_0, irow);
                gtoval::distribute(submat_xz_0, buffer_xz_0, irow);
                gtoval::distribute(submat_yy_0, buffer_yy_0, irow);
                gtoval::distribute(submat_yz_0, buffer_yz_0, irow);
                gtoval::distribute(submat_zz_0, buffer_zz_0, irow);
    
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
get_3rd_order_values_rec_s(const CGtoBlock&            gto_block,
                           const std::vector<double>&  grid_coords_x,
                           const std::vector<double>&  grid_coords_y,
                           const std::vector<double>&  grid_coords_z,
                           const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("3RD_ORDER", nrows, ncols);
        
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
        
        auto submat_0_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x_0 = gto_values.sub_matrix({1, 0});
        auto submat_y_0 = gto_values.sub_matrix({1, 1});
        auto submat_z_0 = gto_values.sub_matrix({1, 2});

        auto submat_xx_0 = gto_values.sub_matrix({2, 0});
        auto submat_xy_0 = gto_values.sub_matrix({2, 1});
        auto submat_xz_0 = gto_values.sub_matrix({2, 2});
        auto submat_yy_0 = gto_values.sub_matrix({2, 3});
        auto submat_yz_0 = gto_values.sub_matrix({2, 4});
        auto submat_zz_0 = gto_values.sub_matrix({2, 5});

        auto submat_xxx_0 = gto_values.sub_matrix({3, 0});
        auto submat_xxy_0 = gto_values.sub_matrix({3, 1});
        auto submat_xxz_0 = gto_values.sub_matrix({3, 2});
        auto submat_xyy_0 = gto_values.sub_matrix({3, 3});
        auto submat_xyz_0 = gto_values.sub_matrix({3, 4});
        auto submat_xzz_0 = gto_values.sub_matrix({3, 5});
        auto submat_yyy_0 = gto_values.sub_matrix({3, 6});
        auto submat_yyz_0 = gto_values.sub_matrix({3, 7});
        auto submat_yzz_0 = gto_values.sub_matrix({3, 8});
        auto submat_zzz_0 = gto_values.sub_matrix({3, 9});

        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer_0_0(ncols);
        
        std::vector<double> buffer_x_0(ncols);
        std::vector<double> buffer_y_0(ncols);
        std::vector<double> buffer_z_0(ncols);

        std::vector<double> buffer_xx_0(ncols);
        std::vector<double> buffer_xy_0(ncols);
        std::vector<double> buffer_xz_0(ncols);
        std::vector<double> buffer_yy_0(ncols);
        std::vector<double> buffer_yz_0(ncols);
        std::vector<double> buffer_zz_0(ncols);

        std::vector<double> buffer_xxx_0(ncols);
        std::vector<double> buffer_xxy_0(ncols);
        std::vector<double> buffer_xxz_0(ncols);
        std::vector<double> buffer_xyy_0(ncols);
        std::vector<double> buffer_xyz_0(ncols);
        std::vector<double> buffer_xzz_0(ncols);
        std::vector<double> buffer_yyy_0(ncols);
        std::vector<double> buffer_yyz_0(ncols);
        std::vector<double> buffer_yzz_0(ncols);
        std::vector<double> buffer_zzz_0(ncols);
        
        auto ptr_buffer_0_0 = buffer_0_0.data();
        
        auto ptr_buffer_x_0 = buffer_x_0.data();
        auto ptr_buffer_y_0 = buffer_y_0.data();
        auto ptr_buffer_z_0 = buffer_z_0.data();

        auto ptr_buffer_xx_0 = buffer_xx_0.data();
        auto ptr_buffer_xy_0 = buffer_xy_0.data();
        auto ptr_buffer_xz_0 = buffer_xz_0.data();
        auto ptr_buffer_yy_0 = buffer_yy_0.data();
        auto ptr_buffer_yz_0 = buffer_yz_0.data();
        auto ptr_buffer_zz_0 = buffer_zz_0.data();

        auto ptr_buffer_xxx_0 = buffer_xxx_0.data();
        auto ptr_buffer_xxy_0 = buffer_xxy_0.data();
        auto ptr_buffer_xxz_0 = buffer_xxz_0.data();
        auto ptr_buffer_xyy_0 = buffer_xyy_0.data();
        auto ptr_buffer_xyz_0 = buffer_xyz_0.data();
        auto ptr_buffer_xzz_0 = buffer_xzz_0.data();
        auto ptr_buffer_yyy_0 = buffer_yyy_0.data();
        auto ptr_buffer_yyz_0 = buffer_yyz_0.data();
        auto ptr_buffer_yzz_0 = buffer_yzz_0.data();
        auto ptr_buffer_zzz_0 = buffer_zzz_0.data();
        
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
                
                std::ranges::fill(buffer_0_0, 0.0);
                
                std::ranges::fill(buffer_x_0, 0.0);
                std::ranges::fill(buffer_y_0, 0.0);
                std::ranges::fill(buffer_z_0, 0.0);

                std::ranges::fill(buffer_xx_0, 0.0);
                std::ranges::fill(buffer_xy_0, 0.0);
                std::ranges::fill(buffer_xz_0, 0.0);
                std::ranges::fill(buffer_yy_0, 0.0);
                std::ranges::fill(buffer_yz_0, 0.0);
                std::ranges::fill(buffer_zz_0, 0.0);

                std::ranges::fill(buffer_xxx_0, 0.0);
                std::ranges::fill(buffer_xxy_0, 0.0);
                std::ranges::fill(buffer_xxz_0, 0.0);
                std::ranges::fill(buffer_xyy_0, 0.0);
                std::ranges::fill(buffer_xyz_0, 0.0);
                std::ranges::fill(buffer_xzz_0, 0.0);
                std::ranges::fill(buffer_yyy_0, 0.0);
                std::ranges::fill(buffer_yyz_0, 0.0);
                std::ranges::fill(buffer_yzz_0, 0.0);
                std::ranges::fill(buffer_zzz_0, 0.0);
              
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

                        const auto fg_0 = -2.0 * fexp;

                        const auto fg_1 = f0_0 * fg_0;

                        const auto fg_2 = fg_1 * fg_0;

                        const auto fg_3 = fg_2 * fg_0;
                        
                        ptr_buffer_0_0[k] += f0_0;
                        
                        ptr_buffer_x_0[k] += gr_x * fg_1;
                        ptr_buffer_y_0[k] += gr_y * fg_1;
                        ptr_buffer_z_0[k] += gr_z * fg_1;

                        ptr_buffer_xx_0[k] += fg_2 * gr_x * gr_x + fg_1;
                        ptr_buffer_xy_0[k] += fg_2 * gr_x * gr_y;
                        ptr_buffer_xz_0[k] += fg_2 * gr_x * gr_z;
                        ptr_buffer_yy_0[k] += fg_2 * gr_y * gr_y + fg_1;
                        ptr_buffer_yz_0[k] += fg_2 * gr_y * gr_z;
                        ptr_buffer_zz_0[k] += fg_2 * gr_z * gr_z + fg_1;

                        ptr_buffer_xxx_0[k] += fg_3 * gr_x * gr_x * gr_x + 3.0 * fg_2 * gr_x;
                        ptr_buffer_xxy_0[k] += fg_3 * gr_x * gr_x * gr_y + fg_2 * gr_y;
                        ptr_buffer_xxz_0[k] += fg_3 * gr_x * gr_x * gr_z + fg_2 * gr_z;
                        ptr_buffer_xyy_0[k] += fg_3 * gr_x * gr_y * gr_y + fg_2 * gr_x;
                        ptr_buffer_xyz_0[k] += fg_3 * gr_x * gr_y * gr_z;
                        ptr_buffer_xzz_0[k] += fg_3 * gr_x * gr_z * gr_z + fg_2 * gr_x;
                        ptr_buffer_yyy_0[k] += fg_3 * gr_y * gr_y * gr_y + 3.0 * fg_2 * gr_y;
                        ptr_buffer_yyz_0[k] += fg_3 * gr_y * gr_y * gr_z + fg_2 * gr_z;
                        ptr_buffer_yzz_0[k] += fg_3 * gr_y * gr_z * gr_z + fg_2 * gr_y;
                        ptr_buffer_zzz_0[k] += fg_3 * gr_z * gr_z * gr_z + 3.0 * fg_2 * gr_z;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0_0, buffer_0_0, irow);
                
                gtoval::distribute(submat_x_0, buffer_x_0, irow);
                gtoval::distribute(submat_y_0, buffer_y_0, irow);
                gtoval::distribute(submat_z_0, buffer_z_0, irow);

                gtoval::distribute(submat_xx_0, buffer_xx_0, irow);
                gtoval::distribute(submat_xy_0, buffer_xy_0, irow);
                gtoval::distribute(submat_xz_0, buffer_xz_0, irow);
                gtoval::distribute(submat_yy_0, buffer_yy_0, irow);
                gtoval::distribute(submat_yz_0, buffer_yz_0, irow);
                gtoval::distribute(submat_zz_0, buffer_zz_0, irow);

                gtoval::distribute(submat_xxx_0, buffer_xxx_0, irow);
                gtoval::distribute(submat_xxy_0, buffer_xxy_0, irow);
                gtoval::distribute(submat_xxz_0, buffer_xxz_0, irow);
                gtoval::distribute(submat_xyy_0, buffer_xyy_0, irow);
                gtoval::distribute(submat_xyz_0, buffer_xyz_0, irow);
                gtoval::distribute(submat_xzz_0, buffer_xzz_0, irow);
                gtoval::distribute(submat_yyy_0, buffer_yyy_0, irow);
                gtoval::distribute(submat_yyz_0, buffer_yyz_0, irow);
                gtoval::distribute(submat_yzz_0, buffer_yzz_0, irow);
                gtoval::distribute(submat_zzz_0, buffer_zzz_0, irow);
    
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
