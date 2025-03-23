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

#include "SerialDensityGridGenerator.hpp"

#include <omp.h>

#include <cstring>
#include <iostream>

#include "MathFunc.hpp"
#include "SerialDenseLinearAlgebra.hpp"

namespace sdengridgen {  // sdengridgen namespace

auto
serialGenerateDensityForLDA(double* rho, const CDenseMatrix& gtoValues, const CDenseMatrix& densityMatrix) -> void
{
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    auto mat_F = sdenblas::serialMultAB(densityMatrix, gtoValues);

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
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    auto mat_F_a = sdenblas::serialMultAB(densityMatrixAlpha, gtoValues);
    auto mat_F_b = sdenblas::serialMultAB(densityMatrixBeta, gtoValues);

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

        auto mat_F = sdenblas::serialMultAB(densityMatrix.getReferenceToDensity(idens), gtoValues);

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

    auto mat_F = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValues);

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

    auto mat_F_a = sdenblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = sdenblas::serialMultAB(symmetricDensityMatrixBeta, gtoValues);

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

auto
serialGenerateDensityGridForGGA(const CDenseMatrix&     gtoValues,
                                const CDenseMatrix&     gtoValuesX,
                                const CDenseMatrix&     gtoValuesY,
                                const CDenseMatrix&     gtoValuesZ,
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

        auto grada = dengrid.alphaDensityGradient(idens);

        auto gradb = dengrid.betaDensityGradient(idens);

        auto gradab = dengrid.mixedDensityGradient(idens);

        auto gradax = dengrid.alphaDensityGradientX(idens);

        auto graday = dengrid.alphaDensityGradientY(idens);

        auto gradaz = dengrid.alphaDensityGradientZ(idens);

        auto gradbx = dengrid.betaDensityGradientX(idens);

        auto gradby = dengrid.betaDensityGradientY(idens);

        auto gradbz = dengrid.betaDensityGradientZ(idens);

        CDenseMatrix symmetricDensityMatrix(densityMatrix.getReferenceToDensity(idens)); 

        symmetricDensityMatrix.symmetrizeAndScale(0.5);

        auto mat_F = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValues);

        auto naos = gtoValues.getNumberOfRows();

        auto F_val = mat_F.values();

        auto chi_val = gtoValues.values();

        auto chi_x_val = gtoValuesX.values();

        auto chi_y_val = gtoValuesY.values();

        auto chi_z_val = gtoValuesZ.values();

        {
            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
                {
                    rhoa[g] += F_val[nu_offset + g] * chi_val[nu_offset + g];

                    gradax[g] += 2.0 * F_val[nu_offset + g] * chi_x_val[nu_offset + g];

                    graday[g] += 2.0 * F_val[nu_offset + g] * chi_y_val[nu_offset + g];

                    gradaz[g] += 2.0 * F_val[nu_offset + g] * chi_z_val[nu_offset + g];
                }
            }

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
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
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValues);

    auto mat_F_x = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesX);
    auto mat_F_y = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesY);
    auto mat_F_z = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesZ);

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
    auto naos = gtoValues.getNumberOfRows();

    auto npoints = gtoValues.getNumberOfColumns();

    CDenseMatrix symmetricDensityMatrixAlpha(densityMatrixAlpha);
    CDenseMatrix symmetricDensityMatrixBeta(densityMatrixBeta);

    symmetricDensityMatrixAlpha.symmetrizeAndScale(0.5);
    symmetricDensityMatrixBeta.symmetrizeAndScale(0.5);

    auto mat_F_a = sdenblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValues);
    auto mat_F_b = sdenblas::serialMultAB(symmetricDensityMatrixBeta, gtoValues);

    auto mat_F_a_x = sdenblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesX);
    auto mat_F_a_y = sdenblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesY);
    auto mat_F_a_z = sdenblas::serialMultAB(symmetricDensityMatrixAlpha, gtoValuesZ);

    auto mat_F_b_x = sdenblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesX);
    auto mat_F_b_y = sdenblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesY);
    auto mat_F_b_z = sdenblas::serialMultAB(symmetricDensityMatrixBeta, gtoValuesZ);

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

auto
serialGenerateDensityGridForMGGA(const CDenseMatrix&     gtoValues,
                                 const CDenseMatrix&     gtoValuesX,
                                 const CDenseMatrix&     gtoValuesY,
                                 const CDenseMatrix&     gtoValuesZ,
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

        CDenseMatrix symmetricDensityMatrix(densityMatrix.getReferenceToDensity(idens)); 

        symmetricDensityMatrix.symmetrizeAndScale(0.5);

        auto mat_F = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValues);

        auto mat_F_x = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesX);
        auto mat_F_y = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesY);
        auto mat_F_z = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValuesZ);

        // eq.(27), JCTC 2021, 17, 1512-1521

        auto naos = gtoValues.getNumberOfRows();

        auto F_val = mat_F.values();

        auto F_x_val = mat_F_x.values();
        auto F_y_val = mat_F_y.values();
        auto F_z_val = mat_F_z.values();

        auto chi_val = gtoValues.values();

        auto chi_x_val = gtoValuesX.values();
        auto chi_y_val = gtoValuesY.values();
        auto chi_z_val = gtoValuesZ.values();

        {
            for (int nu = 0; nu < naos; nu++)
            {
                auto nu_offset = nu * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
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
            for (int g = 0; g < npoints; g++)
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
    }

    return dengrid;
}

}  // namespace sdengridgen
