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

#include "GtoValuesRecS.hpp"

#include <cmath>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
getLdaValuesRecS(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int>(grid_coords_x.size());

    auto gto_values = matfunc::make_matrix("LDA", nrows, ncols);

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

    std::vector<double> buffer(ncols);

    auto ptr_buffer = buffer.data();

    int irow = 0;

    for (int i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto r_x = gto_coords[i].coordinates()[0];

            const auto r_y = gto_coords[i].coordinates()[1];

            const auto r_z = gto_coords[i].coordinates()[2];

            // compute GTO values on grid

            mathfunc::zero(buffer);

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

auto
getGgaValuesRecS(const CGtoBlock&            gto_block,
                 const std::vector<double>&  grid_coords_x,
                 const std::vector<double>&  grid_coords_y,
                 const std::vector<double>&  grid_coords_z,
                 const std::vector<int>& gtos_mask) -> CMatrix
{
    // set up GTO values storage

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int>(grid_coords_x.size());

    auto gto_values = matfunc::make_matrix("GGA", nrows, ncols);

    auto submat_0_0 = gto_values.sub_matrix({0, 0});
    auto submat_x_0 = gto_values.sub_matrix({1, 0});
    auto submat_y_0 = gto_values.sub_matrix({1, 1});
    auto submat_z_0 = gto_values.sub_matrix({1, 2});

    submat_0_0->zero();
    submat_x_0->zero();
    submat_y_0->zero();
    submat_z_0->zero();

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

    std::vector<double> buffer_0_0(ncols);
    std::vector<double> buffer_x_0(ncols);
    std::vector<double> buffer_y_0(ncols);
    std::vector<double> buffer_z_0(ncols);

    auto ptr_buffer_0_0 = buffer_0_0.data();
    auto ptr_buffer_x_0 = buffer_x_0.data();
    auto ptr_buffer_y_0 = buffer_y_0.data();
    auto ptr_buffer_z_0 = buffer_z_0.data();

    int irow = 0;

    for (int i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto r_x = gto_coords[i].coordinates()[0];

            const auto r_y = gto_coords[i].coordinates()[1];

            const auto r_z = gto_coords[i].coordinates()[2];

            // compute GTO values on grid

            mathfunc::zero(buffer_0_0);
            mathfunc::zero(buffer_x_0);
            mathfunc::zero(buffer_y_0);
            mathfunc::zero(buffer_z_0);

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

}  // namespace gtoval
