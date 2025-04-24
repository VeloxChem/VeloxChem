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

#include "PairDensityGridGenerator.hpp"

#include <omp.h>

#include <cstring>

#include "MathFunc.hpp"
#include "SerialDenseLinearAlgebra.hpp"

namespace pairdengridgen {  // pairdengridgen namespace

void
serialGeneratePairDensityForLDA(double*               rho,
                                const CDenseMatrix&   gtoValues,
                                const CDenseMatrix&   densityMatrix,
                                const CDenseMatrix&   activeMOs,
                                const CDenseMatrix&   twoBodyDensityMatrix)
{
    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    auto mat_F = sdenblas::serialMultAB(densityMatrix, gtoValues);

    auto n_active = activeMOs.getNumberOfRows();

    CDenseMatrix MOs_on_grid;

    if (n_active > 0)
    {
        MOs_on_grid = sdenblas::serialMultAB(activeMOs, gtoValues);
    }

    auto n_active2 = n_active * n_active;

    CDenseMatrix mo_pair(n_active2, npoints);

    auto mo_pair_val = mo_pair.values();

    {
        for (int t = 0; t < n_active; t++)
        {
            auto MOt = MOs_on_grid.row(t);

            auto t_offset = t * n_active * npoints;

            for (int u = 0; u < n_active; u++)
            {
                auto MOu = MOs_on_grid.row(u);

                auto tu_offset = t_offset + u * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
                {
                    mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                }
            }
        }
    }

    auto mat_d = sdenblas::serialMultAB(twoBodyDensityMatrix, mo_pair);

    auto naos = gtoValues.getNumberOfRows();

    auto F_val = mat_F.values();

    auto d_val = mat_d.values();

    auto chi_val = gtoValues.values();

    {
        #pragma omp simd 
        for (int g = 0; g < npoints; g++)
        {
            rho[2 * g + 0] = 0.0;
            rho[2 * g + 1] = 0.0;
        }

        // Density

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                rho[2 * g + 0] += F_val[nu_offset + g] * chi_val[nu_offset + g];
            }
        }

        // Pair density

        for (int vw = 0; vw < n_active2; vw++)
        {
            auto vw_offset = vw * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
            {
                rho[2 * g + 1] += d_val[vw_offset + g] * mo_pair_val[vw_offset + g];
            }
        }

        // To prevent numerical issues, enforce that -0.5*rho^2 < pi < 0.5*rho^2

        for (int g = 0; g < npoints; g++)
        {
            auto bound = 0.5 * rho[2 * g + 0] * rho[2 * g + 0];

            rho[2 * g + 1] = std::min(rho[2 * g + 1], bound);

            rho[2 * g + 1] = std::max(rho[2 * g + 1], -bound);
        }
    }
}

void
serialGeneratePairDensityForGGA(double*               rho,
                                double*               rhograd,
                                double*               sigma,
                                const CDenseMatrix&   gtoValues,
                                const CDenseMatrix&   gtoValuesX,
                                const CDenseMatrix&   gtoValuesY,
                                const CDenseMatrix&   gtoValuesZ,
                                const CDenseMatrix&   densityMatrix,
                                const CDenseMatrix&   activeMOs,
                                const CDenseMatrix&   twoBodyDensityMatrix)
{
    auto npoints = gtoValues.getNumberOfColumns();

    // eq.(26), JCTC 2021, 17, 1512-1521

    CDenseMatrix symmetricDensityMatrix(densityMatrix);

    symmetricDensityMatrix.symmetrizeAndScale(0.5);

    auto mat_F = sdenblas::serialMultAB(symmetricDensityMatrix, gtoValues);

    auto n_active = activeMOs.getNumberOfRows();

    CDenseMatrix MOs_on_grid;

    CDenseMatrix MOs_on_gridX;
    CDenseMatrix MOs_on_gridY;
    CDenseMatrix MOs_on_gridZ;

    if (n_active > 0)
    {
        MOs_on_grid = sdenblas::serialMultAB(activeMOs, gtoValues);

        MOs_on_gridX = sdenblas::serialMultAB(activeMOs, gtoValuesX);
        MOs_on_gridY = sdenblas::serialMultAB(activeMOs, gtoValuesY);
        MOs_on_gridZ = sdenblas::serialMultAB(activeMOs, gtoValuesZ);
    }

    auto n_active2 = n_active * n_active;

    CDenseMatrix mo_pair(n_active2, npoints);

    auto mo_pair_val = mo_pair.values();

    {
        for (int t = 0; t < n_active; t++)
        {
            auto MOt = MOs_on_grid.row(t);

            auto t_offset = t * n_active * npoints;

            for (int u = 0; u < n_active; u++)
            {
                auto MOu = MOs_on_grid.row(u);

                auto tu_offset = t_offset + u * npoints;

                #pragma omp simd 
                for (int g = 0; g < npoints; g++)
                {
                    mo_pair_val[tu_offset + g] += MOt[g] * MOu[g];
                }
            }
        }
    }

    auto mat_d = sdenblas::serialMultAB(twoBodyDensityMatrix, mo_pair);

    auto naos = gtoValues.getNumberOfRows();

    auto F_val = mat_F.values();

    auto chi_val = gtoValues.values();

    auto d_val = mat_d.values();

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

        // Density

        for (int nu = 0; nu < naos; nu++)
        {
            auto nu_offset = nu * npoints;

            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
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
            for (int g = 0; g < npoints; g++)
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
                for (int g = 0; g < npoints; g++)
                {
                    rhograd[6 * g + 3] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlX[g];
                    rhograd[6 * g + 4] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlY[g];
                    rhograd[6 * g + 5] += 4.0 * d_val[vw_offset + g] * MOw[g] * MOlZ[g];
                }
            }
        }

        // To prevent numerical issues, enforce that -0.5*rho^2 < pi < 0.5*rho^2

        for (int g = 0; g < npoints; g++)
        {
            auto bound = 0.5 * rho[2 * g + 0] * rho[2 * g + 0];

            rho[2 * g + 1] = std::min(rho[2 * g + 1], bound);

            rho[2 * g + 1] = std::max(rho[2 * g + 1], -bound);
        }

        if (sigma != nullptr)
        {
            #pragma omp simd 
            for (int g = 0; g < npoints; g++)
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
}

}  // namespace pairdengridgen
