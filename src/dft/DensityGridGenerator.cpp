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

#include <cstring>

#include <omp.h>

#include "DenseLinearAlgebra.hpp"

#include <iostream>

namespace dengridgen {  // dengridgen namespace

void
generateDensityForLDA(double*             rho,
                      const int32_t       npoints,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& densityMatrix,
                      CMultiTimer&        timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

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

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_alpha
            rho[2 * g + 0] = 0.0;
        }

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, F_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_alpha
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

        #pragma omp simd aligned(rho : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_beta
            rho[2 * g + 1] = rho[2 * g + 0];
        }
    }

    timer.stop("Density grid rho");
}

void
generateDensityForLDA(double*             rho,
                      const int32_t       npoints,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& densityMatrixAlpha,
                      const CDenseMatrix& densityMatrixBeta,
                      CMultiTimer&        timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F_a = denblas::multAB(densityMatrixAlpha, gtoValues);
    auto mat_F_b = denblas::multAB(densityMatrixBeta, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto chi_val = gtoValues.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;
        }

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, F_a_val, F_b_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }
    }

    timer.stop("Density grid rho");
}

void
generateDensityForGGA(double*             rho,
                      double*             rhograd,
                      double*             sigma,
                      const int32_t       npoints,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& gtoValuesX,
                      const CDenseMatrix& gtoValuesY,
                      const CDenseMatrix& gtoValuesZ,
                      const CDenseMatrix& densityMatrix,
                      CMultiTimer&        timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

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

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho, rhograd : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // rho_a
            rho[2 * g + 0] = 0.0;

            // rho_a_grad
            rhograd[6 * g + 0] = 0.0;
            rhograd[6 * g + 1] = 0.0;
            rhograd[6 * g + 2] = 0.0;
        }

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, rhograd, F_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_a
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                // rho_a_grad
                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        #pragma omp simd aligned(rho, rhograd : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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
            #pragma omp simd aligned(rhograd, sigma : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // simga_aa, sigma_ab, sigma_bb
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

void
generateDensityForGGA(double*             rho,
                      double*             rhograd,
                      double*             sigma,
                      const int32_t       npoints,
                      const CDenseMatrix& gtoValues,
                      const CDenseMatrix& gtoValuesX,
                      const CDenseMatrix& gtoValuesY,
                      const CDenseMatrix& gtoValuesZ,
                      const CDenseMatrix& densityMatrixAlpha,
                      const CDenseMatrix& densityMatrixBeta,
                      CMultiTimer&        timer)
{
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

    auto naos = gtoValues.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

    auto F_a_val = mat_F_a.values();
    auto F_b_val = mat_F_b.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho, rhograd : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, rhograd, F_a_val, F_b_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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
            #pragma omp simd aligned(rhograd, sigma : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

void
generateDensityForMGGA(double*             rho,
                       double*             rhograd,
                       double*             sigma,
                       double*             lapl,
                       double*             tau,
                       const int32_t       npoints,
                       const CDenseMatrix& gtoValues,
                       const CDenseMatrix& gtoValuesX,
                       const CDenseMatrix& gtoValuesY,
                       const CDenseMatrix& gtoValuesZ,
                       const CDenseMatrix& gtoValuesXX,
                       const CDenseMatrix& gtoValuesYY,
                       const CDenseMatrix& gtoValuesZZ,
                       const CDenseMatrix& densityMatrix,
                       CMultiTimer&        timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

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

    auto chi_xx_val = gtoValuesXX.values();
    auto chi_yy_val = gtoValuesYY.values();
    auto chi_zz_val = gtoValuesZZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho, rhograd, lapl, tau : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, rhograd, F_val, chi_val, chi_x_val, chi_y_val, chi_z_val, \
                                     F_x_val, F_y_val, F_z_val, chi_xx_val, chi_yy_val, chi_zz_val, \
                                     lapl, tau : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                // rho_a
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                // rho_a_grad
                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];

                // lapl_a
                lapl[2 * g + 0] += 2.0 * F_val[nu_offset + g] * (chi_xx_val[nu_offset + g] +
                                                                 chi_yy_val[nu_offset + g] +
                                                                 chi_zz_val[nu_offset + g]);

                // tau_a
                tau[2 * g + 0] += 0.5 * (F_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
            }
        }

        #pragma omp simd aligned(rho, rhograd, lapl, tau : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

        #pragma omp simd aligned(rhograd, sigma : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            // simga_aa, sigma_ab, sigma_bb
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

    timer.stop("Density grid rho");
}

void
generateDensityForMGGA(double*             rho,
                       double*             rhograd,
                       double*             sigma,
                       double*             lapl,
                       double*             tau,
                       const int32_t       npoints,
                       const CDenseMatrix& gtoValues,
                       const CDenseMatrix& gtoValuesX,
                       const CDenseMatrix& gtoValuesY,
                       const CDenseMatrix& gtoValuesZ,
                       const CDenseMatrix& gtoValuesXX,
                       const CDenseMatrix& gtoValuesYY,
                       const CDenseMatrix& gtoValuesZZ,
                       const CDenseMatrix& densityMatrixAlpha,
                       const CDenseMatrix& densityMatrixBeta,
                       CMultiTimer&        timer)
{
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

    auto naos = gtoValues.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

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

    auto chi_xx_val = gtoValuesXX.values();
    auto chi_yy_val = gtoValuesYY.values();
    auto chi_zz_val = gtoValuesZZ.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho, rhograd, lapl, tau : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, rhograd, F_a_val, F_b_val, chi_val, \
                                     chi_x_val, chi_y_val, chi_z_val, chi_xx_val, chi_yy_val, chi_zz_val, \
                                     F_a_x_val, F_a_y_val, F_a_z_val, F_b_x_val, F_b_y_val, F_b_z_val, \
                                     lapl, tau : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_a_val[nu_offset + g] * chi_val[nu_offset + g];
                rho[2 * g + 1] += F_b_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_a_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_a_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_a_val[nu_offset + g] * chi_z_val[nu_offset + g];

                rhograd[6 * g + 3] += 2.0 * F_b_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 4] += 2.0 * F_b_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 5] += 2.0 * F_b_val[nu_offset + g] * chi_z_val[nu_offset + g];

                lapl[2 * g + 0] += 2.0 * F_a_val[nu_offset + g] * (chi_xx_val[nu_offset + g] +
                                                                   chi_yy_val[nu_offset + g] +
                                                                   chi_zz_val[nu_offset + g]);
                lapl[2 * g + 1] += 2.0 * F_b_val[nu_offset + g] * (chi_xx_val[nu_offset + g] +
                                                                   chi_yy_val[nu_offset + g] +
                                                                   chi_zz_val[nu_offset + g]);

                tau[2 * g + 0] += 0.5 * (F_a_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_a_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_a_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
                tau[2 * g + 1] += 0.5 * (F_b_x_val[nu_offset + g] * chi_x_val[nu_offset + g] +
                                         F_b_y_val[nu_offset + g] * chi_y_val[nu_offset + g] +
                                         F_b_z_val[nu_offset + g] * chi_z_val[nu_offset + g]);
            }
        }

        #pragma omp simd aligned(rhograd, sigma : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

    timer.stop("Density grid rho");
}

void
generatePairDensityForLDA(double*               rho,
                          const int32_t         npoints,
                          const CDenseMatrix&   gtoValues,
                          const CDenseMatrix&   densityMatrix,
                          const CDenseMatrix&   activeMOs,
                          const CDense4DTensor& twoBodyDensityMatrix,
                          CMultiTimer&          timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

    auto MOs_on_grid = denblas::multAB(activeMOs, gtoValues);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto n_active = activeMOs.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto twoDM = twoBodyDensityMatrix.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;
        }

        // Density

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, F_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

        // Pair density

        for (int32_t i = 0; i < n_active; i++)
        {
            auto MOi = MOs_on_grid.row(i);

            for (int32_t j = 0; j < n_active; j++)
            {
                auto ij = i * n_active + j;

                auto MOj = MOs_on_grid.row(j);

                for (int32_t k = 0; k < n_active; k++)
                {
                    auto ijk = ij * n_active + k;

                    auto MOk = MOs_on_grid.row(k);

                    for (int32_t l = 0; l < n_active; l++)
                    {
                        auto ijkl = ijk * n_active + l;

                        auto twoDM_ijkl = twoDM[ijkl];

                        auto MOl = MOs_on_grid.row(l);

                        #pragma omp simd aligned(rho, MOi, MOj, MOk, MOl : VLX_ALIGN)
                        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                        {
                            rho[2 * g + 1] += twoDM_ijkl * MOi[g] * MOj[g] * MOk[g] * MOl[g];
                        }
                    }
                }
            }
        }
    }

    timer.stop("Density grid rho");
}

void
generatePairDensityForGGA(double*               rho,
                          double*               rhograd,
                          double*               sigma,
                          const int32_t         npoints,
                          const CDenseMatrix&   gtoValues,
                          const CDenseMatrix&   gtoValuesX,
                          const CDenseMatrix&   gtoValuesY,
                          const CDenseMatrix&   gtoValuesZ,
                          const CDenseMatrix&   densityMatrix,
                          const CDenseMatrix&   activeMOs,
                          const CDense4DTensor& twoBodyDensityMatrix,
                          CMultiTimer&          timer)
{
    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = denblas::multAB(symmetricDensityMatrix, gtoValues);

    auto MOs_on_grid = denblas::multAB(activeMOs, gtoValues);

    auto MOs_on_gridX = denblas::multAB(activeMOs, gtoValuesX);
    auto MOs_on_gridY = denblas::multAB(activeMOs, gtoValuesY);
    auto MOs_on_gridZ = denblas::multAB(activeMOs, gtoValuesZ);

    timer.stop("Density grid matmul");

    // eq.(27), JCTC 2021, 17, 1512-1521

    timer.start("Density grid rho");

    auto naos = gtoValues.getNumberOfRows();

    auto n_active = activeMOs.getNumberOfRows();

    auto nthreads = omp_get_max_threads();

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto chi_x_val = gtoValuesX.values();
    auto chi_y_val = gtoValuesY.values();
    auto chi_z_val = gtoValuesZ.values();

    auto twoDM = twoBodyDensityMatrix.values();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        #pragma omp simd aligned(rho, rhograd : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rho, rhograd, F_val, chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                rhograd[6 * g + 0] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];
                rhograd[6 * g + 1] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];
                rhograd[6 * g + 2] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        // Pair density

        for (int32_t i = 0; i < n_active; i++)
        {
            auto MOi = MOs_on_grid.row(i);

            for (int32_t j = 0; j < n_active; j++)
            {
                auto ij = i * n_active + j;

                auto MOj = MOs_on_grid.row(j);

                for (int32_t k = 0; k < n_active; k++)
                {
                    auto ijk = ij * n_active + k;

                    auto MOk = MOs_on_grid.row(k);

                    for (int32_t l = 0; l < n_active; l++)
                    {
                        auto ijkl = ijk * n_active + l;

                        auto twoDM_ijkl = twoDM[ijkl];

                        auto MOl = MOs_on_grid.row(l);

                        auto MOlX = MOs_on_gridX.row(l);
                        auto MOlY = MOs_on_gridY.row(l);
                        auto MOlZ = MOs_on_gridZ.row(l);

                        #pragma omp simd aligned(rho, rhograd, MOi, MOj, MOk, MOl, MOlX, MOlY, MOlZ: VLX_ALIGN)
                        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                        {
                            rho[2 * g + 1] += twoDM_ijkl * MOi[g] * MOj[g] * MOk[g] * MOl[g];

                            rhograd[6 * g + 3] += 4.0 * twoDM_ijkl * MOi[g] * MOj[g] * MOk[g] * MOlX[g];
                            rhograd[6 * g + 4] += 4.0 * twoDM_ijkl * MOi[g] * MOj[g] * MOk[g] * MOlY[g];
                            rhograd[6 * g + 5] += 4.0 * twoDM_ijkl * MOi[g] * MOj[g] * MOk[g] * MOlZ[g];
                        }
                    }
                }
            }
        }

        if (sigma != nullptr)
        {
            #pragma omp simd aligned(rhograd, sigma : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

CDensityGrid
generateDensityGridForLDA(const int32_t       npoints,
                          const CDenseMatrix& gtoValues,
                          const CDenseMatrix& densityMatrix,
                          const xcfun         xcFunType,
                          CMultiTimer&        timer)
{
    CDensityGrid dengrid(npoints, 1, xcFunType, dengrid::ab);

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    auto mat_F = denblas::multAB(densityMatrix, gtoValues);

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

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rhoa, F_val, chi_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

        #pragma omp simd aligned(rhoa, rhob : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
        {
            rhob[g] = rhoa[g];
        }
    }

    timer.stop("Density grid rho");

    return dengrid;
}

CDensityGrid
generateDensityGridForLDA(const int32_t           npoints,
                          const CDenseMatrix&     gtoValues,
                          const CAODensityMatrix& densityMatrix,
                          const xcfun             xcFunType,
                          CMultiTimer&            timer)
{
    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int32_t idens = 0; idens < numdens; idens++)
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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(rhoa, F_val, chi_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];
                }
            }

            #pragma omp simd aligned(rhoa, rhob : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhob[g] = rhoa[g];
            }
        }

        timer.stop("Density grid rho");
    }

    return dengrid;
}

CDensityGrid
generateDensityGridForGGA(const int32_t       npoints,
                          const CDenseMatrix& gtoValues,
                          const CDenseMatrix& gtoValuesX,
                          const CDenseMatrix& gtoValuesY,
                          const CDenseMatrix& gtoValuesZ,
                          const CDenseMatrix& densityMatrix,
                          const xcfun         xcFunType,
                          CMultiTimer&        timer)
{
    CDensityGrid dengrid(npoints, 1, xcFunType, dengrid::ab);

    auto rhoa = dengrid.alphaDensity(0);

    auto rhob = dengrid.betaDensity(0);

    auto grada = dengrid.alphaDensityGradient(0);

    auto gradb = dengrid.betaDensityGradient(0);

    auto gradab = dengrid.mixedDensityGradient(0);

    auto gradax = dengrid.alphaDensityGradientX(0);

    auto graday = dengrid.alphaDensityGradientY(0);

    auto gradaz = dengrid.alphaDensityGradientZ(0);

    auto gradbx = dengrid.betaDensityGradientX(0);

    auto gradby = dengrid.betaDensityGradientY(0);

    auto gradbz = dengrid.betaDensityGradientZ(0);

    // eq.(26), JCTC 2021, 17, 1512-1521

    timer.start("Density grid matmul");

    CDenseMatrix symmetricDensityMatrix(densityMatrix); 

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

        auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

        auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

        for (int32_t nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd aligned(rhoa, gradax, graday, gradaz, F_val, \
                                     chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
            {
                rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];

                graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];

                gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
            }
        }

        #pragma omp simd aligned(rhoa, rhob, gradax, graday, gradaz, gradbx, gradby, gradbz, \
                                 grada, gradb, gradab : VLX_ALIGN)
        for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

    return dengrid;
}

CDensityGrid
generateDensityGridForGGA(const int32_t           npoints,
                          const CDenseMatrix&     gtoValues,
                          const CDenseMatrix&     gtoValuesX,
                          const CDenseMatrix&     gtoValuesY,
                          const CDenseMatrix&     gtoValuesZ,
                          const CAODensityMatrix& densityMatrix,
                          const xcfun             xcFunType,
                          CMultiTimer&            timer)
{
    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    CDensityGrid dengrid(npoints, numdens, xcFunType, dengrid::ab);

    for (int32_t idens = 0; idens < numdens; idens++)
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

            auto grid_batch_size = mpi::batch_size(npoints, thread_id, nthreads);

            auto grid_batch_offset = mpi::batch_offset(npoints, thread_id, nthreads);

            for (int32_t nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd aligned(rhoa, gradax, graday, gradaz, F_val, \
                                         chi_val, chi_x_val, chi_y_val, chi_z_val : VLX_ALIGN)
                for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                    gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];

                    graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];

                    gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
                }
            }

            #pragma omp simd aligned(rhoa, rhob, gradax, graday, gradaz, gradbx, gradby, gradbz, \
                                     grada, gradb, gradab : VLX_ALIGN)
            for (int32_t g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; g++)
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

}  // namespace dengridgen
