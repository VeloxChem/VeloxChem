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

#include "GridScreener.hpp"

namespace gridscreen {  // gridscreen namespace

void
screenDensityGridForLDA(std::vector<int32_t>& gridPointInds,
                        CDensityGrid&         destDensityGrid,
                        const CDensityGrid&   srcDensityGrid,
                        const double          densityThreshold)
{
    auto dest_rhoa = destDensityGrid.alphaDensity(0);

    auto dest_rhob = destDensityGrid.betaDensity(0);

    auto src_rhoa = srcDensityGrid.alphaDensity(0);

    auto src_rhob = srcDensityGrid.betaDensity(0);

    auto npoints = srcDensityGrid.getNumberOfGridPoints();

    int32_t gpcount = 0;

    for (int32_t g = 0; g < npoints; g++)
    {
        if ((std::fabs(src_rhoa[g]) > densityThreshold) &&
            (std::fabs(src_rhob[g]) > densityThreshold))
        {
            gridPointInds[gpcount] = g;

            dest_rhoa[gpcount] = src_rhoa[g];

            dest_rhob[gpcount] = src_rhob[g];

            ++gpcount;
        }
    }

    destDensityGrid.slice(gpcount);
}

void
screenWeights(double*                     screenedWeights,
              const int32_t               gridBlockPosition,
              const double*               weights,
              const std::vector<int32_t>& gridPointInds,
              const int32_t               nScreenedGridPoints)
{
    for (int32_t g = 0; g < nScreenedGridPoints; g++)
    {
        auto g_orig = gridPointInds[g];

        screenedWeights[g] = weights[gridBlockPosition + g_orig];
    }
}

void
screenGtoMatrixForLDA(CDenseMatrix&               screenedGtoValues,
                      const CDenseMatrix&         originalGtoValues,
                      const std::vector<int32_t>& gridPointInds,
                      const int32_t               nScreenedGridPoints)
{
    auto naos = originalGtoValues.getNumberOfRows();

    for (int32_t nu = 0; nu < naos; nu++)
    {
        auto screened_row = screenedGtoValues.row(nu);

        auto original_row = originalGtoValues.row(nu);

        for (int32_t g = 0; g < nScreenedGridPoints; g++)
        {
            auto g_orig = gridPointInds[g];

            screened_row[g] = original_row[g_orig];
        }
    }
}

void
screenDensityGridForGGA(std::vector<int32_t>& gridPointInds,
                        CDensityGrid&         destDensityGrid,
                        const CDensityGrid&   srcDensityGrid,
                        const double          densityThreshold)
{
    auto dest_rhoa = destDensityGrid.alphaDensity(0);

    auto dest_rhob = destDensityGrid.betaDensity(0);

    auto dest_grada = destDensityGrid.alphaDensityGradient(0);

    auto dest_gradb = destDensityGrid.betaDensityGradient(0);

    auto dest_gradab = destDensityGrid.mixedDensityGradient(0);

    auto dest_grada_x = destDensityGrid.alphaDensityGradientX(0);

    auto dest_grada_y = destDensityGrid.alphaDensityGradientY(0);

    auto dest_grada_z = destDensityGrid.alphaDensityGradientZ(0);

    auto dest_gradb_x = destDensityGrid.betaDensityGradientX(0);

    auto dest_gradb_y = destDensityGrid.betaDensityGradientY(0);

    auto dest_gradb_z = destDensityGrid.betaDensityGradientZ(0);

    auto src_rhoa = srcDensityGrid.alphaDensity(0);

    auto src_rhob = srcDensityGrid.betaDensity(0);

    auto src_grada = srcDensityGrid.alphaDensityGradient(0);

    auto src_gradb = srcDensityGrid.betaDensityGradient(0);

    auto src_gradab = srcDensityGrid.mixedDensityGradient(0);

    auto src_grada_x = srcDensityGrid.alphaDensityGradientX(0);

    auto src_grada_y = srcDensityGrid.alphaDensityGradientY(0);

    auto src_grada_z = srcDensityGrid.alphaDensityGradientZ(0);

    auto src_gradb_x = srcDensityGrid.betaDensityGradientX(0);

    auto src_gradb_y = srcDensityGrid.betaDensityGradientY(0);

    auto src_gradb_z = srcDensityGrid.betaDensityGradientZ(0);

    auto npoints = srcDensityGrid.getNumberOfGridPoints();

    int32_t gpcount = 0;

    for (int32_t g = 0; g < npoints; g++)
    {
        if ((std::fabs(src_rhoa[g]) > densityThreshold) &&
            (std::fabs(src_rhob[g]) > densityThreshold) &&
            (std::fabs(src_grada[g]) > densityThreshold) &&
            (std::fabs(src_gradb[g]) > densityThreshold))
        {
            gridPointInds[gpcount] = g;

            dest_rhoa[gpcount] = src_rhoa[g];

            dest_rhob[gpcount] = src_rhob[g];

            dest_grada[gpcount] = src_grada[g];

            dest_gradb[gpcount] = src_gradb[g];

            dest_gradab[gpcount] = src_gradab[g];

            dest_grada_x[gpcount] = src_grada_x[g];

            dest_grada_y[gpcount] = src_grada_y[g];

            dest_grada_z[gpcount] = src_grada_z[g];

            dest_gradb_x[gpcount] = src_gradb_x[g];

            dest_gradb_y[gpcount] = src_gradb_y[g];

            dest_gradb_z[gpcount] = src_gradb_z[g];

            ++gpcount;
        }
    }

    destDensityGrid.slice(gpcount);
}

void
screenGtoMatrixForGGA(CDenseMatrix&               screenedGtoValues,
                      CDenseMatrix&               screenedGtoValuesX,
                      CDenseMatrix&               screenedGtoValuesY,
                      CDenseMatrix&               screenedGtoValuesZ,
                      const CDenseMatrix&         originalGtoValues,
                      const CDenseMatrix&         originalGtoValuesX,
                      const CDenseMatrix&         originalGtoValuesY,
                      const CDenseMatrix&         originalGtoValuesZ,
                      const std::vector<int32_t>& gridPointInds,
                      const int32_t               nScreenedGridPoints)
{
    auto naos = originalGtoValues.getNumberOfRows();

    for (int32_t nu = 0; nu < naos; nu++)
    {
        auto screened_row = screenedGtoValues.row(nu);

        auto screened_x_row = screenedGtoValuesX.row(nu);

        auto screened_y_row = screenedGtoValuesY.row(nu);

        auto screened_z_row = screenedGtoValuesZ.row(nu);

        auto original_row = originalGtoValues.row(nu);

        auto original_x_row = originalGtoValuesX.row(nu);

        auto original_y_row = originalGtoValuesY.row(nu);

        auto original_z_row = originalGtoValuesZ.row(nu);

        for (int32_t g = 0; g < nScreenedGridPoints; g++)
        {
            auto g_orig = gridPointInds[g];

            screened_row[g] = original_row[g_orig];

            screened_x_row[g] = original_x_row[g_orig];

            screened_y_row[g] = original_y_row[g_orig];

            screened_z_row[g] = original_z_row[g_orig];
        }
    }
}

}  // namespace gridscreen
