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

#include "DensityGridGenerator.hpp"

#include <omp.h>

#include <cstring>
#include <iostream>

#include "DenseLinearAlgebra.hpp"
#include "MathFunc.hpp"

namespace dengridgen {  // dengridgen namespace

auto
serialGenerateDensityForLDA(double* rho, const CDenseMatrix& gtoValues, const CDenseMatrix& densityMatrix) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    auto mat_F = denblas::serialMultAB(densityMatrix, gtoValues);

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            // rho_alpha
            rho[2 * g + 0] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                // rho_alpha
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            // rho_beta
            rho[2 * g + 1] = rho[2 * g + 0];
        }
    }
}

auto
serialGenerateDensityForLDA(double*             rho,
                            const CDenseMatrix& gtoValues,
                            const CDenseMatrix& densityMatrixAlpha,
                            const CDenseMatrix& densityMatrixBeta) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    auto mat_F_a = denblas::serialMultAB(densityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::serialMultAB(densityMatrixBeta, gtoValues);

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto chi_val = gtoValues.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }
    }
}

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
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_alpha
            rho[2 * g + 0] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_alpha
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

#pragma omp simd
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_beta
            rho[2 * g + 1] = rho[2 * g + 0];
        }
    }

    timer.stop("Density grid rho");
}

auto
generateDensityForLDA(double*             rho,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& densityMatrixAlpha,
                      const CDenseMatrix& densityMatrixBeta,
                      CMultiTimer&        timer) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F_a = denblas::multAB(densityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::multAB(densityMatrixBeta, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

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

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Density grid rho");
}

auto
serialGenerateDensityGridForLDA(const CDenseMatrix&     gtoValues,
                                const CAODensityMatrix& densityMatrix,
                                const xcfun             xcFunType) -> CDensityGrid
{
    auto npoints = gtoValues.getNumberOfColumns();

    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int idens = 0; idens < numdens; idens++)
    {
        auto rhoa = dengrid.alphaDensity(idens);

        auto rhob = dengrid.betaDensity(idens);

        auto mat_F = denblas::serialMultAB(densityMatrix.getReferenceToDensity(idens), gtoValues);

        auto naos = gtoValues.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_val = gtoValues.values();

        {
            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];
                }
            }

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                rhob[g] = rhoa[g];

                //std::cout << "===rhob===" <<rhob[g]<<std::endl;
            }
        }
    }

    return dengrid;
}

auto
generateDensityGridForLDA(const CDenseMatrix&     gtoValues,
                          const CAODensityMatrix& densityMatrix,
                          const xcfun             xcFunType,
                          CMultiTimer&            timer) -> CDensityGrid
{
    auto npoints = gtoValues.getNumberOfColumns();

    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int idens = 0; idens < numdens; idens++)
    {
        auto rhoa = dengrid.alphaDensity(idens);

        auto rhob = dengrid.betaDensity(idens);

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grid matmul");

        auto mat_F = denblas::multAB(densityMatrix.getReferenceToDensity(idens), gtoValues);

        timer.stop("Density grid matmul");

        // eq.(27), JCTC 2021, 17, 1512-1521

        timer.start("Density grid rho");

        auto naos = gtoValues.getNumberOfRows();

        auto nthreads = omp_get_max_threads();

        auto F_val = mat_F.values();

        auto chi_val = gtoValues.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];
                }
            }

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhob[g] = rhoa[g];

                //std::cout << "===rhob===" <<rhob[g]<<std::endl;
            }
        }

        timer.stop("Density grid rho");
    }

    return dengrid;
}

auto
serialGenerateDensityForGGA(double*             rho,
                            double*             rhograd,
                            double*             sigma,
                            const CDenseMatrix& gtoValues,
                            const CDenseMatrix& gtoValuesX,
                            const CDenseMatrix& gtoValuesY,
                            const CDenseMatrix& gtoValuesZ,
                            const CDenseMatrix& densityMatrix) -> void
{
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::serialMultAB(symmetricDensityMatrix, gtoValues);

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
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
        for (int g = 0; g < npoints; g++)
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
            for (int g = 0; g < npoints; g++)
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
}

auto
serialGenerateDensityForGGA(double*             rho,
                            double*             rhograd,
                            double*             sigma,
                            const CDenseMatrix& gtoValues,
                            const CDenseMatrix& gtoValuesX,
                            const CDenseMatrix& gtoValuesY,
                            const CDenseMatrix& gtoValuesZ,
                            const CDenseMatrix& densityMatrixAlpha,
                            const CDenseMatrix& densityMatrixBeta) -> void
{
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    CDenseMatrix symmetricDensityMatrixAlpha(densityMatrixAlpha);
    CDenseMatrix symmetricDensityMatrixBeta(densityMatrixBeta);

    symmetricDensityMatrixAlpha.symmetrizeAndScale(0.5);
    symmetricDensityMatrixBeta.symmetrizeAndScale(0.5);

    auto mat_F_a = denblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::serialMultAB(symmetricDensityMatrixBeta, gtoValues);

    // eq.(27), JCTC 2021, 17, 1512-1521

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
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

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_a_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_a_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_a_val[nu_offset + g] * chi_z_val[nu_offset + g];
                rhograd[6 * g + 3] += 2.0 * F_b_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 4] += 2.0 * F_b_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 5] += 2.0 * F_b_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        if (sigma != nullptr)
        {
#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
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
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

auto
generateDensityForGGA(double*             rho,
                      double*             rhograd,
                      double*             sigma,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& gtoValuesX,
                      const CDenseMatrix& gtoValuesY,
                      const CDenseMatrix& gtoValuesZ,
                      const CDenseMatrix& densityMatrixAlpha,
                      const CDenseMatrix& densityMatrixBeta,
                      CMultiTimer&        timer) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrixAlpha(densityMatrixAlpha);
    CDenseMatrix symmetricDensityMatrixBeta(densityMatrixBeta);

    symmetricDensityMatrixAlpha.symmetrizeAndScale(0.5);
    symmetricDensityMatrixBeta.symmetrizeAndScale(0.5);

    auto mat_F_a = denblas::multAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::multAB(symmetricDensityMatrixBeta, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

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

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_a_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_a_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_a_val[nu_offset + g] * chi_z_val[nu_offset + g];
                rhograd[6 * g + 3] += 2.0 * F_b_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 4] += 2.0 * F_b_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 5] += 2.0 * F_b_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        if (sigma != nullptr)
        {
#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
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

auto
generateDensityGridForGGA(const CDenseMatrix&     gtoValues,
                          const CDenseMatrix&     gtoValuesX,
                          const CDenseMatrix&     gtoValuesY,
                          const CDenseMatrix&     gtoValuesZ,
                          const CAODensityMatrix& densityMatrix,
                          const xcfun             xcFunType,
                          CMultiTimer&            timer) -> CDensityGrid
{
    auto npoints = gtoValues.getNumberOfColumns();

    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int idens = 0; idens < numdens; idens++)
    {
        auto rhoa = dengrid.alphaDensity(idens);

        auto rhob = dengrid.betaDensity(idens);

        auto grada = dengrid.alphaDensityGradient(idens);

        auto gradb = dengrid.betaDensityGradient(idens);

        auto gradab = dengrid.mixedDensityGradient(idens);

        auto gradax = dengrid.alphaDensityGradientX(idens);

        auto graday = dengrid.alphaDensityGradientY(idens);

        auto gradaz = dengrid.alphaDensityGradientZ(idens);

        auto gradbx = dengrid.betaDensityGradientX(idens);

        auto gradby = dengrid.betaDensityGradientY(idens);

        auto gradbz = dengrid.betaDensityGradientZ(idens);

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grid matmul");

        CDenseMatrix symmetricDensityMatrix(densityMatrix.getReferenceToDensity(idens)); 

        symmetricDensityMatrix.symmetrizeAndScale(0.5);

        auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

        timer.stop("Density grid matmul");

        // eq.(27), JCTC 2021, 17, 1512-1521

        timer.start("Density grid rho");

        auto naos = gtoValues.getNumberOfRows();

        auto nthreads = omp_get_max_threads();

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

            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                    gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];

                    graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];

                    gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
                }
            }

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhob[g] = rhoa[g];

                gradbx[g] = gradax[g];

                gradby[g] = graday[g];

                gradbz[g] = gradaz[g];

                grada[g] = std::sqrt(gradax[g] * gradax[g] + graday[g] * graday[g] + gradaz[g] * gradaz[g]);

                gradb[g] = std::sqrt(gradbx[g] * gradbx[g] + gradby[g] * gradby[g] + gradbz[g] * gradbz[g]);

                gradab[g] = gradax[g] * gradbx[g] + graday[g] * gradby[g] + gradaz[g] * gradbz[g];
            }
        }

        timer.stop("Density grid rho");
    }

    return dengrid;
}

auto
serialGenerateDensityForMGGA(double*             rho,
                             double*             rhograd,
                             double*             sigma,
                             double*             lapl,
                             double*             tau,
                             const CDenseMatrix& gtoValues,
                             const CDenseMatrix& gtoValuesX,
                             const CDenseMatrix& gtoValuesY,
                             const CDenseMatrix& gtoValuesZ,
                             const CDenseMatrix& densityMatrix) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::serialMultAB(symmetricDensityMatrix, gtoValues);

    auto mat_F_x = denblas::serialMultAB(symmetricDensityMatrix, gtoValuesX);
    auto mat_F_y = denblas::serialMultAB(symmetricDensityMatrix, gtoValuesY);
    auto mat_F_z = denblas::serialMultAB(symmetricDensityMatrix, gtoValuesZ);

    auto F_val = mat_F.values();

    auto F_x_val = mat_F_x.values();
    auto F_y_val = mat_F_y.values();
    auto F_z_val = mat_F_z.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;

            // lapl_a
            lapl[2 * g + 0] = 0.0;

            // tau_a
            tau[2 * g + 0] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                // rho_a
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                // rho_a_grad
                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];

                // tau_a

                // clang-format off
                tau[2 * g + 0] += 0.5 * (F_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                // clang-format on
            }
        }

#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            // rho_b
            rho[2 * g + 1] = rho[2 * g + 0];

            // rho_b_grad
            rhograd[6 * g + 3] = rhograd[6 * g + 0];
            rhograd[6 * g + 4] = rhograd[6 * g + 1];
            rhograd[6 * g + 5] = rhograd[6 * g + 2];

            // lapl_b
            lapl[2 * g + 1] = lapl[2 * g + 0];

            // tau_b
            tau[2 * g + 1] = tau[2 * g + 0];
        }

        if (sigma != nullptr)
        {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
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
}

auto
serialGenerateDensityForMGGA(double*             rho,
                             double*             rhograd,
                             double*             sigma,
                             double*             lapl,
                             double*             tau,
                             const CDenseMatrix& gtoValues,
                             const CDenseMatrix& gtoValuesX,
                             const CDenseMatrix& gtoValuesY,
                             const CDenseMatrix& gtoValuesZ,
                             const CDenseMatrix& densityMatrixAlpha,
                             const CDenseMatrix& densityMatrixBeta) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    CDenseMatrix symmetricDensityMatrixAlpha(densityMatrixAlpha);
    CDenseMatrix symmetricDensityMatrixBeta(densityMatrixBeta);

    symmetricDensityMatrixAlpha.symmetrizeAndScale(0.5);
    symmetricDensityMatrixBeta.symmetrizeAndScale(0.5);

    auto mat_F_a = denblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::serialMultAB(symmetricDensityMatrixBeta, gtoValues);

    auto mat_F_a_x = denblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesX);
    auto mat_F_a_y = denblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesY);
    auto mat_F_a_z = denblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesZ);

    auto mat_F_b_x = denblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesX);
    auto mat_F_b_y = denblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesY);
    auto mat_F_b_z = denblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesZ);

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto F_a_x_val = mat_F_a_x.values();
    auto F_a_y_val = mat_F_a_y.values();
    auto F_a_z_val = mat_F_a_z.values();

    auto F_b_x_val = mat_F_b_x.values();
    auto F_b_y_val = mat_F_b_y.values();
    auto F_b_z_val = mat_F_b_z.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;

            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
            rhograd[6 * g + 3] = 0.0;
            rhograd[6 * g + 4] = 0.0;
            rhograd[6 * g + 5] = 0.0;

            lapl[2 * g + 0] = 0.0;
            lapl[2 * g + 1] = 0.0;

            tau[2 * g + 0] = 0.0;
            tau[2 * g + 1] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = 0; g < npoints; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_a_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_a_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_a_val[nu_offset + g] * chi_z_val[nu_offset + g];

                rhograd[6 * g + 3] += 2.0 * F_b_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 4] += 2.0 * F_b_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 5] += 2.0 * F_b_val[nu_offset + g] * chi_z_val[nu_offset + g];

                // TODO implement Laplacian dependence

                // clang-format off
                tau[2 * g + 0] += 0.5 * (F_a_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_a_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_a_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);

                tau[2 * g + 1] += 0.5 * (F_b_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_b_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_b_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                // clang-format on
            }
        }

        if (sigma != nullptr)
        {
#pragma omp simd
        for (int g = 0; g < npoints; g++)
        {
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
}

void
generateDensityForMGGA(double*             rho,
                       double*             rhograd,
                       double*             sigma,
                       double*             lapl,
                       double*             tau,
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

    auto mat_F_x = denblas::multAB(symmetricDensityMatrix, gtoValuesX);
    auto mat_F_y = denblas::multAB(symmetricDensityMatrix, gtoValuesY);
    auto mat_F_z = denblas::multAB(symmetricDensityMatrix, gtoValuesZ);

    timer.stop("Density grid matmul");

    timer.start("Density grid rho");

    auto F_val = mat_F.values();

    auto F_x_val = mat_F_x.values();
    auto F_y_val = mat_F_y.values();
    auto F_z_val = mat_F_z.values();

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
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;

            // lapl_a
            lapl[2 * g + 0] = 0.0;

            // tau_a
            tau[2 * g + 0] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_a
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                // rho_a_grad
                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];

                // tau_a

                // clang-format off
                tau[2 * g + 0] += 0.5 * (F_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                // clang-format on
            }
        }

#pragma omp simd
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_b
            rho[2 * g + 1] = rho[2 * g + 0];

            // rho_b_grad
            rhograd[6 * g + 3] = rhograd[6 * g + 0];
            rhograd[6 * g + 4] = rhograd[6 * g + 1];
            rhograd[6 * g + 5] = rhograd[6 * g + 2];

            // lapl_b
            lapl[2 * g + 1] = lapl[2 * g + 0];

            // tau_b
            tau[2 * g + 1] = tau[2 * g + 0];
        }

#pragma omp simd
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

    timer.stop("Density grid rho");
}

auto
generateDensityForMGGA(double*             rho,
                       double*             rhograd,
                       double*             sigma,
                       double*             lapl,
                       double*             tau,
                       const CDenseMatrix& gtoValues,
                       const CDenseMatrix& gtoValuesX,
                       const CDenseMatrix& gtoValuesY,
                       const CDenseMatrix& gtoValuesZ,
                       const CDenseMatrix& densityMatrixAlpha,
                       const CDenseMatrix& densityMatrixBeta,
                       CMultiTimer&        timer) -> void
{
    auto nthreads = omp_get_max_threads();

    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrixAlpha(densityMatrixAlpha);
    CDenseMatrix symmetricDensityMatrixBeta(densityMatrixBeta);

    symmetricDensityMatrixAlpha.symmetrizeAndScale(0.5);
    symmetricDensityMatrixBeta.symmetrizeAndScale(0.5);

    auto mat_F_a = denblas::multAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::multAB(symmetricDensityMatrixBeta, gtoValues);

    auto mat_F_a_x = denblas::multAB(symmetricDensityMatrixAlpha, gtoValuesX);
    auto mat_F_a_y = denblas::multAB(symmetricDensityMatrixAlpha, gtoValuesY);
    auto mat_F_a_z = denblas::multAB(symmetricDensityMatrixAlpha, gtoValuesZ);

    auto mat_F_b_x = denblas::multAB(symmetricDensityMatrixBeta, gtoValuesX);
    auto mat_F_b_y = denblas::multAB(symmetricDensityMatrixBeta, gtoValuesY);
    auto mat_F_b_z = denblas::multAB(symmetricDensityMatrixBeta, gtoValuesZ);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto F_a_x_val = mat_F_a_x.values();
    auto F_a_y_val = mat_F_a_y.values();
    auto F_a_z_val = mat_F_a_z.values();

    auto F_b_x_val = mat_F_b_x.values();
    auto F_b_y_val = mat_F_b_y.values();
    auto F_b_z_val = mat_F_b_z.values();

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

            lapl[2 * g + 0] = 0.0;
            lapl[2 * g + 1] = 0.0;

            tau[2 * g + 0] = 0.0;
            tau[2 * g + 1] = 0.0;
        }

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

#pragma omp simd
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_a_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_a_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_a_val[nu_offset + g] * chi_z_val[nu_offset + g];

                rhograd[6 * g + 3] += 2.0 * F_b_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 4] += 2.0 * F_b_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 5] += 2.0 * F_b_val[nu_offset + g] * chi_z_val[nu_offset + g];

                // TODO implement Laplacian dependence

                // clang-format off
                tau[2 * g + 0] += 0.5 * (F_a_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_a_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_a_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);

                tau[2 * g + 1] += 0.5 * (F_b_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_b_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_b_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                // clang-format on
            }
        }

#pragma omp simd
        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
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

    timer.stop("Density grid rho");
}

auto
generateDensityGridForMGGA(const CDenseMatrix&     gtoValues,
                           const CDenseMatrix&     gtoValuesX,
                           const CDenseMatrix&     gtoValuesY,
                           const CDenseMatrix&     gtoValuesZ,
                           const CAODensityMatrix& densityMatrix,
                           const xcfun             xcFunType,
                           CMultiTimer&            timer) -> CDensityGrid
{
    auto npoints = gtoValues.getNumberOfColumns();

    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int idens = 0; idens < numdens; idens++)
    {
        auto rhoa = dengrid.alphaDensity(idens);
        auto rhob = dengrid.betaDensity(idens);

        auto taua = dengrid.alphaDensitytau(idens);
        auto taub = dengrid.betaDensitytau(idens);

        auto lapla = dengrid.alphaDensitylapl(idens);
        auto laplb = dengrid.betaDensitylapl(idens);

        auto grada = dengrid.alphaDensityGradient(idens);
        auto gradb = dengrid.betaDensityGradient(idens);
        auto gradab = dengrid.mixedDensityGradient(idens);

        auto gradax = dengrid.alphaDensityGradientX(idens);
        auto graday = dengrid.alphaDensityGradientY(idens);
        auto gradaz = dengrid.alphaDensityGradientZ(idens);

        auto gradbx = dengrid.betaDensityGradientX(idens);
        auto gradby = dengrid.betaDensityGradientY(idens);
        auto gradbz = dengrid.betaDensityGradientZ(idens);

        // eq.(26), JCTC 2021, 17, 1512-1521

        timer.start("Density grid matmul");

        CDenseMatrix symmetricDensityMatrix(densityMatrix.getReferenceToDensity(idens)); 

        symmetricDensityMatrix.symmetrizeAndScale(0.5);

        auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

        auto mat_F_x = denblas::multAB(symmetricDensityMatrix, gtoValuesX);
        auto mat_F_y = denblas::multAB(symmetricDensityMatrix, gtoValuesY);
        auto mat_F_z = denblas::multAB(symmetricDensityMatrix, gtoValuesZ);

        timer.stop("Density grid matmul");

        // eq.(27), JCTC 2021, 17, 1512-1521

        timer.start("Density grid rho");

        auto naos = gtoValues.getNumberOfRows();

        auto nthreads = omp_get_max_threads();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_val = gtoValues.values();

        auto chi_x_val = gtoValuesX.values();
        auto chi_y_val = gtoValuesY.values();
        auto chi_z_val = gtoValuesZ.values();

        #pragma omp parallel
        {
            auto thread_id = omp_get_thread_num();

            auto grid_batch_size = mathfunc::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mathfunc::batch_offset(npoints, thread_id, nthreads);

            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                    gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                    graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                    gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];

                    // TODO implement Laplacian dependence

                    // tau_a
                    taua[g] += 0.5 * (F_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                      F_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                      F_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                }
            }

            #pragma omp simd 
            for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhob[g] = rhoa[g];
                laplb[g] = lapla[g];
                taub[g] = taua[g]; 

                gradbx[g] = gradax[g];
                gradby[g] = graday[g];
                gradbz[g] = gradaz[g];

                grada[g] = std::sqrt(gradax[g] * gradax[g] + graday[g] * graday[g] + gradaz[g] * gradaz[g]);
                gradb[g] = std::sqrt(gradbx[g] * gradbx[g] + gradby[g] * gradby[g] + gradbz[g] * gradbz[g]);
                gradab[g] = gradax[g] * gradbx[g] + graday[g] * gradby[g] + gradaz[g] * gradbz[g];
            }
        }

        timer.stop("Density grid rho");
    }

    return dengrid;
}

}  // namespace dengridgen
