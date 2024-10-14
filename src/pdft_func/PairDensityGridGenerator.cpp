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

#include "PairDensityGridGenerator.hpp"

#include <omp.h>

#include <cstring>

#include "DenseLinearAlgebra.hpp"
#include "MathFunc.hpp"

namespace pairdengridgen {  // pairdengridgen namespace

void
generatePairDensityForLDA(double*               rho,
                          const CDenseMatrix&   gtoValues,
                          const CDenseMatrix&   densityMatrix,
                          const CDenseMatrix&   activeMOs,
                          const CDenseMatrix&   twoBodyDensityMatrix,
                          CMultiTimer&          timer)
{
    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

    auto n_active = activeMOs.getNumberOfRows();

    CDenseMatrix MOs_on_grid;
    if (n_active > 0)
    {
        MOs_on_grid = denblas::multAB(activeMOs, gtoValues);
    }

    timer.stop("Density grid matmul");

    timer.start("Density mo pair");

    auto n_active2 = n_active * n_active;

    auto nthreads = omp_get_max_threads();

    CDenseMatrix mo_pair(n_active2, npoints);

    auto mo_pair_val = mo_pair.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int t = 0; t < n_active; t++)
        {
            auto MOt = MOs_on_grid.row(t);

            auto t_offset = t * n_active * npoints;

            for (int u = 0; u < n_active; u++)
            {
                auto MOu = MOs_on_grid.row(u);

                auto tu_offset = t_offset + u * npoints;
                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                }
            }
        }
    }

    timer.stop("Density mo pair");

    timer.start("Density grid pi matmul");

    auto mat_d = denblas::multAB(twoBodyDensityMatrix, mo_pair);

    timer.stop("Density grid pi matmul");
    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto F_val = mat_F.values();

    auto d_val = mat_d.values();

    auto chi_val = gtoValues.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd 
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;
        }

        // Density

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

        // Pair density

        for (int vw = 0; vw < n_active2; vw++)
        {
            auto vw_offset = vw * npoints;

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 1] += d_val[vw_offset + g] * mo_pair_val[vw_offset + g];
            }
        }

        // To prevent numerical issues, enforce that -0.5*rho^2 < pi < 0.5*rho^2

        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            auto bound = 0.5 * rho[2 * g + 0] * rho[2 * g + 0];

            rho[2 * g + 1] = std::min(rho[2 * g + 1], bound);

            rho[2 * g + 1] = std::max(rho[2 * g + 1], -bound);
        }
    }


    timer.stop("Density grid rho");
}

void
generatePairDensityForGGA(double*               rho,
                          double*               rhograd,
                          double*               sigma,
                          const CDenseMatrix&   gtoValues,
                          const CDenseMatrix&   gtoValuesX,
                          const CDenseMatrix&   gtoValuesY,
                          const CDenseMatrix&   gtoValuesZ,
                          const CDenseMatrix&   densityMatrix,
                          const CDenseMatrix&   activeMOs,
                          const CDenseMatrix&   twoBodyDensityMatrix,
                          CMultiTimer&          timer)
{
    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

    auto n_active = activeMOs.getNumberOfRows();

    CDenseMatrix MOs_on_grid;

    CDenseMatrix MOs_on_gridX;
    CDenseMatrix MOs_on_gridY;
    CDenseMatrix MOs_on_gridZ;

    if (n_active > 0)
    {
        MOs_on_grid = denblas::multAB(activeMOs, gtoValues);

        MOs_on_gridX = denblas::multAB(activeMOs, gtoValuesX);
        MOs_on_gridY = denblas::multAB(activeMOs, gtoValuesY);
        MOs_on_gridZ = denblas::multAB(activeMOs, gtoValuesZ);
    }

    timer.stop("Density grid matmul");

    timer.start("Density mo pair");

    auto n_active2 = n_active * n_active;

    auto nthreads = omp_get_max_threads();

    CDenseMatrix mo_pair(n_active2, npoints);

    auto mo_pair_val = mo_pair.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        for (int t = 0; t < n_active; t++)
        {
            auto MOt = MOs_on_grid.row(t);

            auto t_offset = t * n_active * npoints;

            for (int u = 0; u < n_active; u++)
            {
                auto MOu = MOs_on_grid.row(u);

                auto tu_offset = t_offset + u * npoints;
                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                }
            }
        }
    }

    timer.stop("Density mo pair");

    timer.start("Density grid pi matmul");

    auto mat_d = denblas::multAB(twoBodyDensityMatrix, mo_pair);

    timer.stop("Density grid pi matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto d_val = mat_d.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd 
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;

            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
            rhograd[6 * g + 3] = 0.0;
            rhograd[6 * g + 4] = 0.0;
            rhograd[6 * g + 5] = 0.0;
        }

        // Density

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        // Pair density

        for (int vw = 0; vw < n_active2; vw++)
        {
            auto vw_offset = vw * npoints;

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 1] += d_val[vw_offset + g] * mo_pair_val[vw_offset + g];
            }
        }

        for (int v = 0; v < n_active; v++)
        {
            auto MOlX = MOs_on_gridX.row(v);
            auto MOlY = MOs_on_gridY.row(v);
            auto MOlZ = MOs_on_gridZ.row(v);

            for (int w = 0; w < n_active; w++)
            {
                auto vw_offset = (v * n_active + w)* npoints;

                auto MOw = MOs_on_grid.row(w);

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhograd[6 * g + 3] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlX[g];
                    rhograd[6 * g + 4] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlY[g];
                    rhograd[6 * g + 5] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlZ[g];
                }
            }
        }

        // To prevent numerical issues, enforce that -0.5*rho^2 < pi < 0.5*rho^2

        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            auto bound = 0.5 * rho[2 * g + 0] * rho[2 * g + 0];

            rho[2 * g + 1] = std::min(rho[2 * g + 1], bound);

            rho[2 * g + 1] = std::max(rho[2 * g + 1], -bound);
        }

        if (sigma != nullptr)
        {
            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                sigma[3 * g + 0] = rhograd[6 * g + 0] * rhograd[6 * g + 0] +
                                   rhograd[6 * g + 1] * rhograd[6 * g + 1] +
                                   rhograd[6 * g + 2] * rhograd[6 * g + 2];

                sigma[3 * g + 1] = rhograd[6 * g + 0] * rhograd[6 * g + 3] +
                                   rhograd[6 * g + 1] * rhograd[6 * g + 4] +
                                   rhograd[6 * g + 2] * rhograd[6 * g + 5];

                sigma[3 * g + 2] = rhograd[6 * g + 3] * rhograd[6 * g + 3] +
                                   rhograd[6 * g + 4] * rhograd[6 * g + 4] +
                                   rhograd[6 * g + 5] * rhograd[6 * g + 5];
            }
        }
    }

    timer.stop("Density grid rho");
}

}  // namespace pairdengridgen
