//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include "DensityGridGenerator.hpp"

#include <omp.h>

#include <cstring>
#include <iostream>

#include "DenseLinearAlgebra.hpp"
#include "MathFunc.hpp"

namespace dengridgen {  // dengridgen namespace

auto
generateDensityForLDA(double* rho, const CDenseMatrix& gtoValues, const CDenseMatrix& densityMatrix, CMultiTimer& timer) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

    timer.stop("Density grid matmul");

    timer.start("Density grid rho");

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

#pragma omp simd
        for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_alpha
            rho[2 * g + 0] = 0.0;
        }

        for (int64_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_alpha
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

#pragma omp simd
        for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_beta
            rho[2 * g + 1] = rho[2 * g + 0];
        }
    }

    timer.stop("Density grid rho");
}

void
generateDensityForGGA(double*             rho,
                      double*             rhograd,
                      double*             sigma,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& gtoValuesX,
                      const CDenseMatrix& gtoValuesY,
                      const CDenseMatrix& gtoValuesZ,
                      const CDenseMatrix& densityMatrix,
                      CMultiTimer&        timer)
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

    timer.stop("Density grid matmul");

    timer.start("Density grid rho");

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

#pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

#pragma omp simd
        for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
        }

        for (int64_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_a
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                // rho_a_grad
                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

#pragma omp simd
        for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_b
            rho[2 * g + 1] = rho[2 * g + 0];

            // rho_b_grad
            rhograd[6 * g + 3] = rhograd[6 * g + 0];
            rhograd[6 * g + 4] = rhograd[6 * g + 1];
            rhograd[6 * g + 5] = rhograd[6 * g + 2];
        }

        if (sigma != nullptr)
        {
#pragma omp simd
            for (int64_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // simga_aa, sigma_ab, sigma_bb
                // clang-format off
                sigma[3 * g + 0] = rhograd[6 * g + 0] * rhograd[6 * g + 0] +
                                   rhograd[6 * g + 1] * rhograd[6 * g + 1] +
                                   rhograd[6 * g + 2] * rhograd[6 * g + 2];

                sigma[3 * g + 1] = rhograd[6 * g + 0] * rhograd[6 * g + 3] +
                                   rhograd[6 * g + 1] * rhograd[6 * g + 4] +
                                   rhograd[6 * g + 2] * rhograd[6 * g + 5];

                sigma[3 * g + 2] = rhograd[6 * g + 3] * rhograd[6 * g + 3] +
                                   rhograd[6 * g + 4] * rhograd[6 * g + 4] +
                                   rhograd[6 * g + 5] * rhograd[6 * g + 5];
                // clang-format on
            }
        }
    }

    timer.stop("Density grid rho");
}

}  // namespace dengridgen
