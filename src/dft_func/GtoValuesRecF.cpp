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
getLdaValuesRecF(const CGtoBlock&            gto_block,
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

    const auto nrows = mathfunc::countSignificantElements(gtos_mask);

    const auto ncols = static_cast<int>(grid_coords_x.size());

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

            const auto r_x = gto_coords[i].coordinates()[0];

            const auto r_y = gto_coords[i].coordinates()[1];

            const auto r_z = gto_coords[i].coordinates()[2];

            // compute GTO values on grid

            mathfunc::zero(buffer_xxx);

            mathfunc::zero(buffer_xxy);

            mathfunc::zero(buffer_xxz);

            mathfunc::zero(buffer_xyy);

            mathfunc::zero(buffer_xyz);

            mathfunc::zero(buffer_xzz);

            mathfunc::zero(buffer_yyy);

            mathfunc::zero(buffer_yyz);

            mathfunc::zero(buffer_yzz);

            mathfunc::zero(buffer_zzz);

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

}  // namespace gtoval
