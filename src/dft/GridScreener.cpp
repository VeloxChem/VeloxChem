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
screenVxcFockForLDA(double*       rho,
                    double*       exc,
                    double*       vrho,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            vrho[2 * g + 0] = 0.0;
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            vrho[2 * g + 1] = 0.0;
        }

        // rho
        if (std::fabs(rho[2 * g + 0] + rho[2 * g + 1]) <= densityThreshold)
        {
            exc[g] = 0.0;
        }
    }
}

void
screenVxcFockForGGA(double*       rho,
                    double*       sigma,
                    double*       exc,
                    double*       vrho,
                    double*       vsigma,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    double densityThresholdSquared = densityThreshold * densityThreshold;

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 0]) <= densityThresholdSquared))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 2]) <= densityThresholdSquared))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;
        }

        // rho
        if (std::fabs(rho[2 * g + 0] + rho[2 * g + 1]) <= densityThreshold)
        {
            exc[g] = 0.0;
        }
    }
}

void
screenVxcFockForMGGA(double*       rho,
                     double*       sigma,
                     double*       exc,
                     double*       vrho,
                     double*       vsigma,
                     double*       vlapl,
                     double*       vtau,
                     const int32_t npoints,
                     const double  densityThreshold)
{
    double densityThresholdSquared = densityThreshold * densityThreshold;

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 0]) <= densityThresholdSquared))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;

            vlapl[2 * g + 0] = 0.0;

            vtau[2 * g + 0] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 2]) <= densityThresholdSquared))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            vlapl[2 * g + 1] = 0.0;

            vtau[2 * g + 1] = 0.0;
        }

        // rho
        if (std::fabs(rho[2 * g + 0] + rho[2 * g + 1]) <= densityThreshold)
        {
            exc[g] = 0.0;
        }
    }
}

void
screenVxcFockForPLDA(double*       rho,
                     double*       exc,
                     double*       vrho,
                     const int32_t npoints,
                     const double  densityThreshold)
{
    for (int32_t g = 0; g < npoints; g++)
    {
        // rho
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            exc[g] = 0.0;
        }
    }
}

void
screenFxcFockForLDA(double*       rho,
                    double*       v2rho2,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;
        }
    }
}

void
screenFxcFockForGGA(double*       rho,
                    double*       sigma,
                    double*       vrho,
                    double*       vsigma,
                    double*       v2rho2,
                    double*       v2rhosigma,
                    double*       v2sigma2,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    double densityThresholdSquared = densityThreshold * densityThreshold;

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 0]) <= densityThresholdSquared))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;

            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;

            v2rhosigma[6 * g + 0] = 0.0;
            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;

            v2sigma2[6 * g + 0] = 0.0;
            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 2]) <= densityThresholdSquared))
        {
            vrho[2 * g + 1] = 0.0;

            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;
            v2rhosigma[6 * g + 5] = 0.0;

            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;
            v2sigma2[6 * g + 5] = 0.0;
        }
    }
}

void
screenKxcFockForLDA(double*       rho,
                    double*       v2rho2,
                    double*       v3rho3,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;

            v3rho3[4 * g + 0] = 0.0;
            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;

            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;
            v3rho3[4 * g + 3] = 0.0;
        }
    }
}

void
screenKxcFockForGGA(double*       rho,
                    double*       sigma,
                    double*       vrho,
                    double*       vsigma,
                    double*       v2rho2,
                    double*       v2rhosigma,
                    double*       v2sigma2,
                    double*       v3rho3,
                    double*       v3rho2sigma,
                    double*       v3rhosigma2,
                    double*       v3sigma3,
                    const int32_t npoints,
                    const double  densityThreshold)
{
    double densityThresholdSquared = densityThreshold * densityThreshold;

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 0]) <= densityThresholdSquared))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;

            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;

            v2rhosigma[6 * g + 0] = 0.0;
            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;

            v2sigma2[6 * g + 0] = 0.0;
            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;

            v3rho3[4 * g + 0] = 0.0;
            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;

            v3rho2sigma[9 * g + 0] = 0.0;
            v3rho2sigma[9 * g + 1] = 0.0;
            v3rho2sigma[9 * g + 2] = 0.0;
            v3rho2sigma[9 * g + 3] = 0.0;
            v3rho2sigma[9 * g + 4] = 0.0;
            v3rho2sigma[9 * g + 5] = 0.0;
            v3rho2sigma[9 * g + 6] = 0.0;
            v3rho2sigma[9 * g + 7] = 0.0;

            v3rhosigma2[12 * g + 0] = 0.0;
            v3rhosigma2[12 * g + 1] = 0.0;
            v3rhosigma2[12 * g + 2] = 0.0;
            v3rhosigma2[12 * g + 3] = 0.0;
            v3rhosigma2[12 * g + 4] = 0.0;
            v3rhosigma2[12 * g + 5] = 0.0;
            v3rhosigma2[12 * g + 6] = 0.0;
            v3rhosigma2[12 * g + 7] = 0.0;
            v3rhosigma2[12 * g + 8] = 0.0;
            v3rhosigma2[12 * g + 9] = 0.0;
            v3rhosigma2[12 * g + 10] = 0.0;

            v3sigma3[10 * g + 0] = 0.0;
            v3sigma3[10 * g + 1] = 0.0;
            v3sigma3[10 * g + 2] = 0.0;
            v3sigma3[10 * g + 3] = 0.0;
            v3sigma3[10 * g + 4] = 0.0;
            v3sigma3[10 * g + 5] = 0.0;
            v3sigma3[10 * g + 6] = 0.0;
            v3sigma3[10 * g + 7] = 0.0;
            v3sigma3[10 * g + 8] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) ||
            (std::fabs(sigma[3 * g + 2]) <= densityThresholdSquared))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;

            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;
            v2rhosigma[6 * g + 5] = 0.0;

            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;
            v2sigma2[6 * g + 5] = 0.0;

            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;
            v3rho3[4 * g + 3] = 0.0;

            v3rho2sigma[9 * g + 1] = 0.0;
            v3rho2sigma[9 * g + 2] = 0.0;
            v3rho2sigma[9 * g + 3] = 0.0;
            v3rho2sigma[9 * g + 4] = 0.0;
            v3rho2sigma[9 * g + 5] = 0.0;
            v3rho2sigma[9 * g + 6] = 0.0;
            v3rho2sigma[9 * g + 7] = 0.0;
            v3rho2sigma[9 * g + 8] = 0.0;

            v3rhosigma2[12 * g + 1] = 0.0;
            v3rhosigma2[12 * g + 2] = 0.0;
            v3rhosigma2[12 * g + 3] = 0.0;
            v3rhosigma2[12 * g + 4] = 0.0;
            v3rhosigma2[12 * g + 5] = 0.0;
            v3rhosigma2[12 * g + 6] = 0.0;
            v3rhosigma2[12 * g + 7] = 0.0;
            v3rhosigma2[12 * g + 8] = 0.0;
            v3rhosigma2[12 * g + 9] = 0.0;
            v3rhosigma2[12 * g + 10] = 0.0;
            v3rhosigma2[12 * g + 11] = 0.0;

            v3sigma3[10 * g + 1] = 0.0;
            v3sigma3[10 * g + 2] = 0.0;
            v3sigma3[10 * g + 3] = 0.0;
            v3sigma3[10 * g + 4] = 0.0;
            v3sigma3[10 * g + 5] = 0.0;
            v3sigma3[10 * g + 6] = 0.0;
            v3sigma3[10 * g + 7] = 0.0;
            v3sigma3[10 * g + 8] = 0.0;
            v3sigma3[10 * g + 9] = 0.0;
        }
    }
}

void
copyWeights(double*       screenedWeights,
            const int32_t gridBlockPosition,
            const double* weights,
            const int32_t nGridPoints)
{
    for (int32_t g = 0; g < nGridPoints; g++)
    {
        screenedWeights[g] = weights[gridBlockPosition + g];
    }
}

}  // namespace gridscreen
