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

}  // namespace gtoval
