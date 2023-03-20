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

double
getDensityScreeningThreshold()
{
    return 1.0e-15;
}

double
getSigmaScreeningThreshold(const double densityThreshold)
{
    return 1.0e-20;  // std::pow(densityThreshold, 4.0 / 3.0)
}

double
getTauScreeningThreshold()
{
    return 1.0e-20;
}

void
screenExcVxcForLDA(const int32_t npoints, const double* rho, double* exc, double* vrho)
{
    double densityThreshold = getDensityScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho
        if (((std::fabs(rho[2 * g + 0]) <= densityThreshold) && (std::fabs(rho[2 * g + 1]) <= densityThreshold)))
        {
            exc[g] = 0.0;
        }

        // rho_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold))
        {
            vrho[2 * g + 0] = 0.0;
        }

        // rho_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold))
        {
            vrho[2 * g + 1] = 0.0;
        }
    }
}

void
screenVxcForLDA(const int32_t npoints, const double* rho, double* vrho)
{
    double densityThreshold = getDensityScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold))
        {
            vrho[2 * g + 0] = 0.0;
        }

        // rho_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold))
        {
            vrho[2 * g + 1] = 0.0;
        }
    }
}

void
screenExcVxcForGGA(const int32_t npoints, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho and sigma
        if (((std::fabs(rho[2 * g + 0]) <= densityThreshold) && (std::fabs(rho[2 * g + 1]) <= densityThreshold)) ||
            ((std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) && (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold)))
        {
            exc[g] = 0.0;
        }

        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;
        }
    }
}

void
screenVxcForGGA(const int32_t npoints, const double* rho, const double* sigma, double* vrho, double* vsigma)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;
        }
    }
}

void
screenExcVxcForMGGA(const int32_t npoints,
                    const double* rho,
                    const double* sigma,
                    const double* lapl,
                    const double* tau,
                    double*       exc,
                    double*       vrho,
                    double*       vsigma,
                    double*       vlapl,
                    double*       vtau)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho, tau and sigma
        if (((std::fabs(rho[2 * g + 0]) <= densityThreshold) && (std::fabs(rho[2 * g + 1]) <= densityThreshold)) ||
            ((std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) && (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold)) ||
            ((std::fabs(tau[2 * g + 0]) <= tauThreshold) && (std::fabs(tau[2 * g + 1]) <= tauThreshold)))
        {
            exc[g] = 0.0;
        }

        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;

            vlapl[2 * g + 0] = 0.0;

            vtau[2 * g + 0] = 0.0;
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            vlapl[2 * g + 1] = 0.0;

            vtau[2 * g + 1] = 0.0;
        }
    }
}

void
screenVxcForMGGA(const int32_t npoints,
                 const double* rho,
                 const double* sigma,
                 const double* lapl,
                 const double* tau,
                 double*       vrho,
                 double*       vsigma,
                 double*       vlapl,
                 double*       vtau)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            vrho[2 * g + 0] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;

            vlapl[2 * g + 0] = 0.0;

            vtau[2 * g + 0] = 0.0;
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            vlapl[2 * g + 1] = 0.0;

            vtau[2 * g + 1] = 0.0;
        }
    }
}

void
screenVxcForPLDA(const int32_t npoints, const double* rho, double* exc, double* vrho)
{
    double densityThreshold = getDensityScreeningThreshold();

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
screenVxcForPGGA(const int32_t npoints, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double densityThreshold = getDensityScreeningThreshold();

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
screenFxcForLDA(const int32_t npoints, const double* rho, double* v2rho2)
{
    double densityThreshold = getDensityScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold))
        {
            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;
        }

        // rho_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold))
        {
            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;
        }
    }
}

void
screenFxcForGGA(const int32_t npoints, const double* rho, const double* sigma, double* v2rho2, double* v2rhosigma, double* v2sigma2)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
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
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
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
        }
    }
}

void
screenFxcForMGGA(const int32_t npoints,
                 const double* rho,
                 const double* sigma,
                 const double* lapl,
                 const double* tau,
                 double*       v2rho2,
                 double*       v2rhosigma,
                 double*       v2rholapl,
                 double*       v2rhotau,
                 double*       v2sigma2,
                 double*       v2sigmalapl,
                 double*       v2sigmatau,
                 double*       v2lapl2,
                 double*       v2lapltau,
                 double*       v2tau2)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            v2rho2[3 * g + 0] = 0.0;
            v2rho2[3 * g + 1] = 0.0;

            v2rhosigma[6 * g + 0] = 0.0;
            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;

            v2rholapl[4 * g + 0] = 0.0;
            v2rholapl[4 * g + 1] = 0.0;
            v2rholapl[4 * g + 2] = 0.0;

            v2rhotau[4 * g + 0] = 0.0;
            v2rhotau[4 * g + 1] = 0.0;
            v2rhotau[4 * g + 2] = 0.0;

            v2sigma2[6 * g + 0] = 0.0;
            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;

            v2sigmalapl[6 * g + 0] = 0.0;
            v2sigmalapl[6 * g + 1] = 0.0;
            v2sigmalapl[6 * g + 2] = 0.0;
            v2sigmalapl[6 * g + 3] = 0.0;
            v2sigmalapl[6 * g + 4] = 0.0;

            v2sigmatau[6 * g + 0] = 0.0;
            v2sigmatau[6 * g + 1] = 0.0;
            v2sigmatau[6 * g + 2] = 0.0;
            v2sigmatau[6 * g + 3] = 0.0;
            v2sigmatau[6 * g + 4] = 0.0;

            v2lapl2[3 * g + 0] = 0.0;
            v2lapl2[3 * g + 1] = 0.0;

            v2lapltau[4 * g + 0] = 0.0;
            v2lapltau[4 * g + 1] = 0.0;
            v2lapltau[4 * g + 2] = 0.0;

            v2tau2[3 * g + 0] = 0.0;
            v2tau2[3 * g + 1] = 0.0;
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            v2rho2[3 * g + 1] = 0.0;
            v2rho2[3 * g + 2] = 0.0;

            v2rhosigma[6 * g + 1] = 0.0;
            v2rhosigma[6 * g + 2] = 0.0;
            v2rhosigma[6 * g + 3] = 0.0;
            v2rhosigma[6 * g + 4] = 0.0;
            v2rhosigma[6 * g + 5] = 0.0;

            v2rholapl[4 * g + 1] = 0.0;
            v2rholapl[4 * g + 2] = 0.0;
            v2rholapl[4 * g + 3] = 0.0;

            v2rhotau[4 * g + 1] = 0.0;
            v2rhotau[4 * g + 2] = 0.0;
            v2rhotau[4 * g + 3] = 0.0;

            v2sigma2[6 * g + 1] = 0.0;
            v2sigma2[6 * g + 2] = 0.0;
            v2sigma2[6 * g + 3] = 0.0;
            v2sigma2[6 * g + 4] = 0.0;
            v2sigma2[6 * g + 5] = 0.0;

            v2sigmalapl[6 * g + 1] = 0.0;
            v2sigmalapl[6 * g + 2] = 0.0;
            v2sigmalapl[6 * g + 3] = 0.0;
            v2sigmalapl[6 * g + 4] = 0.0;
            v2sigmalapl[6 * g + 5] = 0.0;

            v2sigmatau[6 * g + 1] = 0.0;
            v2sigmatau[6 * g + 2] = 0.0;
            v2sigmatau[6 * g + 3] = 0.0;
            v2sigmatau[6 * g + 4] = 0.0;
            v2sigmatau[6 * g + 5] = 0.0;

            v2lapl2[3 * g + 1] = 0.0;
            v2lapl2[3 * g + 2] = 0.0;

            v2lapltau[4 * g + 1] = 0.0;
            v2lapltau[4 * g + 2] = 0.0;
            v2lapltau[4 * g + 3] = 0.0;

            v2tau2[3 * g + 1] = 0.0;
            v2tau2[3 * g + 2] = 0.0;
        }
    }
}

void
screenKxcForLDA(const int32_t npoints, const double* rho, double* v3rho3)
{
    double densityThreshold = getDensityScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold))
        {
            v3rho3[4 * g + 0] = 0.0;
            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;
        }

        // rho_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold))
        {
            v3rho3[4 * g + 1] = 0.0;
            v3rho3[4 * g + 2] = 0.0;
            v3rho3[4 * g + 3] = 0.0;
        }
    }
}

void
screenKxcForGGA(const int32_t npoints,
                const double* rho,
                const double* sigma,
                double*       v3rho3,
                double*       v3rho2sigma,
                double*       v3rhosigma2,
                double*       v3sigma3)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
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

            v3rhosigma2[12 * g + 0]  = 0.0;
            v3rhosigma2[12 * g + 1]  = 0.0;
            v3rhosigma2[12 * g + 2]  = 0.0;
            v3rhosigma2[12 * g + 3]  = 0.0;
            v3rhosigma2[12 * g + 4]  = 0.0;
            v3rhosigma2[12 * g + 5]  = 0.0;
            v3rhosigma2[12 * g + 6]  = 0.0;
            v3rhosigma2[12 * g + 7]  = 0.0;
            v3rhosigma2[12 * g + 8]  = 0.0;
            v3rhosigma2[12 * g + 9]  = 0.0;
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
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
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

            v3rhosigma2[12 * g + 1]  = 0.0;
            v3rhosigma2[12 * g + 2]  = 0.0;
            v3rhosigma2[12 * g + 3]  = 0.0;
            v3rhosigma2[12 * g + 4]  = 0.0;
            v3rhosigma2[12 * g + 5]  = 0.0;
            v3rhosigma2[12 * g + 6]  = 0.0;
            v3rhosigma2[12 * g + 7]  = 0.0;
            v3rhosigma2[12 * g + 8]  = 0.0;
            v3rhosigma2[12 * g + 9]  = 0.0;
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
screenKxcForMGGA(const int32_t npoints,
                 const double* rho,
                 const double* sigma,
                 const double* lapl,
                 const double* tau,
                 double*       v3rho3,
                 double*       v3rho2sigma,
                 double*       v3rho2lapl,
                 double*       v3rho2tau,
                 double*       v3rhosigma2,
                 double*       v3rhosigmalapl,
                 double*       v3rhosigmatau,
                 double*       v3rholapl2,
                 double*       v3rholapltau,
                 double*       v3rhotau2,
                 double*       v3sigma3,
                 double*       v3sigma2lapl,
                 double*       v3sigma2tau,
                 double*       v3sigmalapl2,
                 double*       v3sigmalapltau,
                 double*       v3sigmatau2,
                 double*       v3lapl3,
                 double*       v3lapl2tau,
                 double*       v3lapltau2,
                 double*       v3tau3)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
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

            v3rho2lapl[6 * g + 0] = 0.0;
            v3rho2lapl[6 * g + 1] = 0.0;
            v3rho2lapl[6 * g + 2] = 0.0;
            v3rho2lapl[6 * g + 3] = 0.0;
            v3rho2lapl[6 * g + 4] = 0.0;

            v3rho2tau[6 * g + 0] = 0.0;
            v3rho2tau[6 * g + 1] = 0.0;
            v3rho2tau[6 * g + 2] = 0.0;
            v3rho2tau[6 * g + 3] = 0.0;
            v3rho2tau[6 * g + 4] = 0.0;

            v3rhosigma2[12 * g + 0]  = 0.0;
            v3rhosigma2[12 * g + 1]  = 0.0;
            v3rhosigma2[12 * g + 2]  = 0.0;
            v3rhosigma2[12 * g + 3]  = 0.0;
            v3rhosigma2[12 * g + 4]  = 0.0;
            v3rhosigma2[12 * g + 5]  = 0.0;
            v3rhosigma2[12 * g + 6]  = 0.0;
            v3rhosigma2[12 * g + 7]  = 0.0;
            v3rhosigma2[12 * g + 8]  = 0.0;
            v3rhosigma2[12 * g + 9]  = 0.0;
            v3rhosigma2[12 * g + 10] = 0.0;

            v3rhosigmalapl[12 * g + 0]  = 0.0;
            v3rhosigmalapl[12 * g + 1]  = 0.0;
            v3rhosigmalapl[12 * g + 2]  = 0.0;
            v3rhosigmalapl[12 * g + 3]  = 0.0;
            v3rhosigmalapl[12 * g + 4]  = 0.0;
            v3rhosigmalapl[12 * g + 5]  = 0.0;
            v3rhosigmalapl[12 * g + 6]  = 0.0;
            v3rhosigmalapl[12 * g + 7]  = 0.0;
            v3rhosigmalapl[12 * g + 8]  = 0.0;
            v3rhosigmalapl[12 * g + 9]  = 0.0;
            v3rhosigmalapl[12 * g + 10] = 0.0;

            v3rhosigmatau[12 * g + 0]  = 0.0;
            v3rhosigmatau[12 * g + 1]  = 0.0;
            v3rhosigmatau[12 * g + 2]  = 0.0;
            v3rhosigmatau[12 * g + 3]  = 0.0;
            v3rhosigmatau[12 * g + 4]  = 0.0;
            v3rhosigmatau[12 * g + 5]  = 0.0;
            v3rhosigmatau[12 * g + 6]  = 0.0;
            v3rhosigmatau[12 * g + 7]  = 0.0;
            v3rhosigmatau[12 * g + 8]  = 0.0;
            v3rhosigmatau[12 * g + 9]  = 0.0;
            v3rhosigmatau[12 * g + 10] = 0.0;

            v3rholapl2[6 * g + 0] = 0.0;
            v3rholapl2[6 * g + 1] = 0.0;
            v3rholapl2[6 * g + 2] = 0.0;
            v3rholapl2[6 * g + 3] = 0.0;
            v3rholapl2[6 * g + 4] = 0.0;

            v3rholapltau[8 * g + 0] = 0.0;
            v3rholapltau[8 * g + 1] = 0.0;
            v3rholapltau[8 * g + 2] = 0.0;
            v3rholapltau[8 * g + 3] = 0.0;
            v3rholapltau[8 * g + 4] = 0.0;
            v3rholapltau[8 * g + 5] = 0.0;
            v3rholapltau[8 * g + 6] = 0.0;

            v3rhotau2[6 * g + 0] = 0.0;
            v3rhotau2[6 * g + 1] = 0.0;
            v3rhotau2[6 * g + 2] = 0.0;
            v3rhotau2[6 * g + 3] = 0.0;
            v3rhotau2[6 * g + 4] = 0.0;

            v3sigma3[10 * g + 0] = 0.0;
            v3sigma3[10 * g + 1] = 0.0;
            v3sigma3[10 * g + 2] = 0.0;
            v3sigma3[10 * g + 3] = 0.0;
            v3sigma3[10 * g + 4] = 0.0;
            v3sigma3[10 * g + 5] = 0.0;
            v3sigma3[10 * g + 6] = 0.0;
            v3sigma3[10 * g + 7] = 0.0;
            v3sigma3[10 * g + 8] = 0.0;

            v3sigma2lapl[12 * g + 0]  = 0.0;
            v3sigma2lapl[12 * g + 1]  = 0.0;
            v3sigma2lapl[12 * g + 2]  = 0.0;
            v3sigma2lapl[12 * g + 3]  = 0.0;
            v3sigma2lapl[12 * g + 4]  = 0.0;
            v3sigma2lapl[12 * g + 5]  = 0.0;
            v3sigma2lapl[12 * g + 6]  = 0.0;
            v3sigma2lapl[12 * g + 7]  = 0.0;
            v3sigma2lapl[12 * g + 8]  = 0.0;
            v3sigma2lapl[12 * g + 9]  = 0.0;
            v3sigma2lapl[12 * g + 10] = 0.0;

            v3sigma2tau[12 * g + 0]  = 0.0;
            v3sigma2tau[12 * g + 1]  = 0.0;
            v3sigma2tau[12 * g + 2]  = 0.0;
            v3sigma2tau[12 * g + 3]  = 0.0;
            v3sigma2tau[12 * g + 4]  = 0.0;
            v3sigma2tau[12 * g + 5]  = 0.0;
            v3sigma2tau[12 * g + 6]  = 0.0;
            v3sigma2tau[12 * g + 7]  = 0.0;
            v3sigma2tau[12 * g + 8]  = 0.0;
            v3sigma2tau[12 * g + 9]  = 0.0;
            v3sigma2tau[12 * g + 10] = 0.0;

            v3sigmalapl2[9 * g + 0] = 0.0;
            v3sigmalapl2[9 * g + 1] = 0.0;
            v3sigmalapl2[9 * g + 2] = 0.0;
            v3sigmalapl2[9 * g + 3] = 0.0;
            v3sigmalapl2[9 * g + 4] = 0.0;
            v3sigmalapl2[9 * g + 5] = 0.0;
            v3sigmalapl2[9 * g + 6] = 0.0;
            v3sigmalapl2[9 * g + 7] = 0.0;

            v3sigmalapltau[12 * g + 0]  = 0.0;
            v3sigmalapltau[12 * g + 1]  = 0.0;
            v3sigmalapltau[12 * g + 2]  = 0.0;
            v3sigmalapltau[12 * g + 3]  = 0.0;
            v3sigmalapltau[12 * g + 4]  = 0.0;
            v3sigmalapltau[12 * g + 5]  = 0.0;
            v3sigmalapltau[12 * g + 6]  = 0.0;
            v3sigmalapltau[12 * g + 7]  = 0.0;
            v3sigmalapltau[12 * g + 8]  = 0.0;
            v3sigmalapltau[12 * g + 9]  = 0.0;
            v3sigmalapltau[12 * g + 10] = 0.0;

            v3sigmatau2[9 * g + 0] = 0.0;
            v3sigmatau2[9 * g + 1] = 0.0;
            v3sigmatau2[9 * g + 2] = 0.0;
            v3sigmatau2[9 * g + 3] = 0.0;
            v3sigmatau2[9 * g + 4] = 0.0;
            v3sigmatau2[9 * g + 5] = 0.0;
            v3sigmatau2[9 * g + 6] = 0.0;
            v3sigmatau2[9 * g + 7] = 0.0;

            v3lapl3[4 * g + 0] = 0.0;
            v3lapl3[4 * g + 1] = 0.0;
            v3lapl3[4 * g + 2] = 0.0;

            v3lapl2tau[6 * g + 0] = 0.0;
            v3lapl2tau[6 * g + 1] = 0.0;
            v3lapl2tau[6 * g + 2] = 0.0;
            v3lapl2tau[6 * g + 3] = 0.0;
            v3lapl2tau[6 * g + 4] = 0.0;

            v3lapltau2[6 * g + 0] = 0.0;
            v3lapltau2[6 * g + 1] = 0.0;
            v3lapltau2[6 * g + 2] = 0.0;
            v3lapltau2[6 * g + 3] = 0.0;
            v3lapltau2[6 * g + 4] = 0.0;

            v3tau3[4 * g + 0] = 0.0;
            v3tau3[4 * g + 1] = 0.0;
            v3tau3[4 * g + 2] = 0.0;
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
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

            v3rho2lapl[6 * g + 1] = 0.0;
            v3rho2lapl[6 * g + 2] = 0.0;
            v3rho2lapl[6 * g + 3] = 0.0;
            v3rho2lapl[6 * g + 4] = 0.0;
            v3rho2lapl[6 * g + 5] = 0.0;

            v3rho2tau[6 * g + 1] = 0.0;
            v3rho2tau[6 * g + 2] = 0.0;
            v3rho2tau[6 * g + 3] = 0.0;
            v3rho2tau[6 * g + 4] = 0.0;
            v3rho2tau[6 * g + 5] = 0.0;

            v3rhosigma2[12 * g + 1]  = 0.0;
            v3rhosigma2[12 * g + 2]  = 0.0;
            v3rhosigma2[12 * g + 3]  = 0.0;
            v3rhosigma2[12 * g + 4]  = 0.0;
            v3rhosigma2[12 * g + 5]  = 0.0;
            v3rhosigma2[12 * g + 6]  = 0.0;
            v3rhosigma2[12 * g + 7]  = 0.0;
            v3rhosigma2[12 * g + 8]  = 0.0;
            v3rhosigma2[12 * g + 9]  = 0.0;
            v3rhosigma2[12 * g + 10] = 0.0;
            v3rhosigma2[12 * g + 11] = 0.0;

            v3rhosigmalapl[12 * g + 1]  = 0.0;
            v3rhosigmalapl[12 * g + 2]  = 0.0;
            v3rhosigmalapl[12 * g + 3]  = 0.0;
            v3rhosigmalapl[12 * g + 4]  = 0.0;
            v3rhosigmalapl[12 * g + 5]  = 0.0;
            v3rhosigmalapl[12 * g + 6]  = 0.0;
            v3rhosigmalapl[12 * g + 7]  = 0.0;
            v3rhosigmalapl[12 * g + 8]  = 0.0;
            v3rhosigmalapl[12 * g + 9]  = 0.0;
            v3rhosigmalapl[12 * g + 10] = 0.0;
            v3rhosigmalapl[12 * g + 11] = 0.0;

            v3rhosigmatau[12 * g + 1]  = 0.0;
            v3rhosigmatau[12 * g + 2]  = 0.0;
            v3rhosigmatau[12 * g + 3]  = 0.0;
            v3rhosigmatau[12 * g + 4]  = 0.0;
            v3rhosigmatau[12 * g + 5]  = 0.0;
            v3rhosigmatau[12 * g + 6]  = 0.0;
            v3rhosigmatau[12 * g + 7]  = 0.0;
            v3rhosigmatau[12 * g + 8]  = 0.0;
            v3rhosigmatau[12 * g + 9]  = 0.0;
            v3rhosigmatau[12 * g + 10] = 0.0;
            v3rhosigmatau[12 * g + 11] = 0.0;

            v3rholapl2[6 * g + 1] = 0.0;
            v3rholapl2[6 * g + 2] = 0.0;
            v3rholapl2[6 * g + 3] = 0.0;
            v3rholapl2[6 * g + 4] = 0.0;
            v3rholapl2[6 * g + 5] = 0.0;

            v3rholapltau[8 * g + 1] = 0.0;
            v3rholapltau[8 * g + 2] = 0.0;
            v3rholapltau[8 * g + 3] = 0.0;
            v3rholapltau[8 * g + 4] = 0.0;
            v3rholapltau[8 * g + 5] = 0.0;
            v3rholapltau[8 * g + 6] = 0.0;
            v3rholapltau[8 * g + 7] = 0.0;

            v3rhotau2[6 * g + 1] = 0.0;
            v3rhotau2[6 * g + 2] = 0.0;
            v3rhotau2[6 * g + 3] = 0.0;
            v3rhotau2[6 * g + 4] = 0.0;
            v3rhotau2[6 * g + 5] = 0.0;

            v3sigma3[10 * g + 1] = 0.0;
            v3sigma3[10 * g + 2] = 0.0;
            v3sigma3[10 * g + 3] = 0.0;
            v3sigma3[10 * g + 4] = 0.0;
            v3sigma3[10 * g + 5] = 0.0;
            v3sigma3[10 * g + 6] = 0.0;
            v3sigma3[10 * g + 7] = 0.0;
            v3sigma3[10 * g + 8] = 0.0;
            v3sigma3[10 * g + 9] = 0.0;

            v3sigma2lapl[12 * g + 1]  = 0.0;
            v3sigma2lapl[12 * g + 2]  = 0.0;
            v3sigma2lapl[12 * g + 3]  = 0.0;
            v3sigma2lapl[12 * g + 4]  = 0.0;
            v3sigma2lapl[12 * g + 5]  = 0.0;
            v3sigma2lapl[12 * g + 6]  = 0.0;
            v3sigma2lapl[12 * g + 7]  = 0.0;
            v3sigma2lapl[12 * g + 8]  = 0.0;
            v3sigma2lapl[12 * g + 9]  = 0.0;
            v3sigma2lapl[12 * g + 10] = 0.0;
            v3sigma2lapl[12 * g + 11] = 0.0;

            v3sigma2tau[12 * g + 1]  = 0.0;
            v3sigma2tau[12 * g + 2]  = 0.0;
            v3sigma2tau[12 * g + 3]  = 0.0;
            v3sigma2tau[12 * g + 4]  = 0.0;
            v3sigma2tau[12 * g + 5]  = 0.0;
            v3sigma2tau[12 * g + 6]  = 0.0;
            v3sigma2tau[12 * g + 7]  = 0.0;
            v3sigma2tau[12 * g + 8]  = 0.0;
            v3sigma2tau[12 * g + 9]  = 0.0;
            v3sigma2tau[12 * g + 10] = 0.0;
            v3sigma2tau[12 * g + 11] = 0.0;

            v3sigmalapl2[9 * g + 1] = 0.0;
            v3sigmalapl2[9 * g + 2] = 0.0;
            v3sigmalapl2[9 * g + 3] = 0.0;
            v3sigmalapl2[9 * g + 4] = 0.0;
            v3sigmalapl2[9 * g + 5] = 0.0;
            v3sigmalapl2[9 * g + 6] = 0.0;
            v3sigmalapl2[9 * g + 7] = 0.0;
            v3sigmalapl2[9 * g + 8] = 0.0;

            v3sigmalapltau[12 * g + 1]  = 0.0;
            v3sigmalapltau[12 * g + 2]  = 0.0;
            v3sigmalapltau[12 * g + 3]  = 0.0;
            v3sigmalapltau[12 * g + 4]  = 0.0;
            v3sigmalapltau[12 * g + 5]  = 0.0;
            v3sigmalapltau[12 * g + 6]  = 0.0;
            v3sigmalapltau[12 * g + 7]  = 0.0;
            v3sigmalapltau[12 * g + 8]  = 0.0;
            v3sigmalapltau[12 * g + 9]  = 0.0;
            v3sigmalapltau[12 * g + 10] = 0.0;
            v3sigmalapltau[12 * g + 11] = 0.0;

            v3sigmatau2[9 * g + 1] = 0.0;
            v3sigmatau2[9 * g + 2] = 0.0;
            v3sigmatau2[9 * g + 3] = 0.0;
            v3sigmatau2[9 * g + 4] = 0.0;
            v3sigmatau2[9 * g + 5] = 0.0;
            v3sigmatau2[9 * g + 6] = 0.0;
            v3sigmatau2[9 * g + 7] = 0.0;
            v3sigmatau2[9 * g + 8] = 0.0;

            v3lapl3[4 * g + 1] = 0.0;
            v3lapl3[4 * g + 2] = 0.0;
            v3lapl3[4 * g + 3] = 0.0;

            v3lapl2tau[6 * g + 1] = 0.0;
            v3lapl2tau[6 * g + 2] = 0.0;
            v3lapl2tau[6 * g + 3] = 0.0;
            v3lapl2tau[6 * g + 4] = 0.0;
            v3lapl2tau[6 * g + 5] = 0.0;

            v3lapltau2[6 * g + 1] = 0.0;
            v3lapltau2[6 * g + 2] = 0.0;
            v3lapltau2[6 * g + 3] = 0.0;
            v3lapltau2[6 * g + 4] = 0.0;
            v3lapltau2[6 * g + 5] = 0.0;

            v3tau3[4 * g + 1] = 0.0;
            v3tau3[4 * g + 2] = 0.0;
            v3tau3[4 * g + 3] = 0.0;
        }
    }
}

void
screenLxcForLDA(const int32_t npoints, const double* rho, double* v4rho4)
{
    double densityThreshold = getDensityScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold))
        {
            v4rho4[5 * g + 0] = 0.0;
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;
        }

        // rho_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold))
        {
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;
            v4rho4[5 * g + 4] = 0.0;
        }
    }
}

void
screenLxcForGGA(const int32_t npoints,
                const double* rho,
                const double* sigma,
                double*       v4rho4,
                double*       v4rho3sigma,
                double*       v4rho2sigma2,
                double*       v4rhosigma3,
                double*       v4sigma4)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            v4rho4[5 * g + 0] = 0.0;
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;

            v4rho3sigma[12 * g + 0]  = 0.0;
            v4rho3sigma[12 * g + 1]  = 0.0;
            v4rho3sigma[12 * g + 2]  = 0.0;
            v4rho3sigma[12 * g + 3]  = 0.0;
            v4rho3sigma[12 * g + 4]  = 0.0;
            v4rho3sigma[12 * g + 5]  = 0.0;
            v4rho3sigma[12 * g + 6]  = 0.0;
            v4rho3sigma[12 * g + 7]  = 0.0;
            v4rho3sigma[12 * g + 8]  = 0.0;
            v4rho3sigma[12 * g + 9]  = 0.0;
            v4rho3sigma[12 * g + 10] = 0.0;

            v4rho2sigma2[18 * g + 0]  = 0.0;
            v4rho2sigma2[18 * g + 1]  = 0.0;
            v4rho2sigma2[18 * g + 2]  = 0.0;
            v4rho2sigma2[18 * g + 3]  = 0.0;
            v4rho2sigma2[18 * g + 4]  = 0.0;
            v4rho2sigma2[18 * g + 5]  = 0.0;
            v4rho2sigma2[18 * g + 6]  = 0.0;
            v4rho2sigma2[18 * g + 7]  = 0.0;
            v4rho2sigma2[18 * g + 8]  = 0.0;
            v4rho2sigma2[18 * g + 9]  = 0.0;
            v4rho2sigma2[18 * g + 10] = 0.0;
            v4rho2sigma2[18 * g + 11] = 0.0;
            v4rho2sigma2[18 * g + 12] = 0.0;
            v4rho2sigma2[18 * g + 13] = 0.0;
            v4rho2sigma2[18 * g + 14] = 0.0;
            v4rho2sigma2[18 * g + 15] = 0.0;
            v4rho2sigma2[18 * g + 16] = 0.0;

            v4rhosigma3[20 * g + 0]  = 0.0;
            v4rhosigma3[20 * g + 1]  = 0.0;
            v4rhosigma3[20 * g + 2]  = 0.0;
            v4rhosigma3[20 * g + 3]  = 0.0;
            v4rhosigma3[20 * g + 4]  = 0.0;
            v4rhosigma3[20 * g + 5]  = 0.0;
            v4rhosigma3[20 * g + 6]  = 0.0;
            v4rhosigma3[20 * g + 7]  = 0.0;
            v4rhosigma3[20 * g + 8]  = 0.0;
            v4rhosigma3[20 * g + 9]  = 0.0;
            v4rhosigma3[20 * g + 10] = 0.0;
            v4rhosigma3[20 * g + 11] = 0.0;
            v4rhosigma3[20 * g + 12] = 0.0;
            v4rhosigma3[20 * g + 13] = 0.0;
            v4rhosigma3[20 * g + 14] = 0.0;
            v4rhosigma3[20 * g + 15] = 0.0;
            v4rhosigma3[20 * g + 16] = 0.0;
            v4rhosigma3[20 * g + 17] = 0.0;
            v4rhosigma3[20 * g + 18] = 0.0;

            v4sigma4[15 * g + 0]  = 0.0;
            v4sigma4[15 * g + 1]  = 0.0;
            v4sigma4[15 * g + 2]  = 0.0;
            v4sigma4[15 * g + 3]  = 0.0;
            v4sigma4[15 * g + 4]  = 0.0;
            v4sigma4[15 * g + 5]  = 0.0;
            v4sigma4[15 * g + 6]  = 0.0;
            v4sigma4[15 * g + 7]  = 0.0;
            v4sigma4[15 * g + 8]  = 0.0;
            v4sigma4[15 * g + 9]  = 0.0;
            v4sigma4[15 * g + 10] = 0.0;
            v4sigma4[15 * g + 11] = 0.0;
            v4sigma4[15 * g + 12] = 0.0;
            v4sigma4[15 * g + 13] = 0.0;
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;
            v4rho4[5 * g + 4] = 0.0;

            v4rho3sigma[12 * g + 1]  = 0.0;
            v4rho3sigma[12 * g + 2]  = 0.0;
            v4rho3sigma[12 * g + 3]  = 0.0;
            v4rho3sigma[12 * g + 4]  = 0.0;
            v4rho3sigma[12 * g + 5]  = 0.0;
            v4rho3sigma[12 * g + 6]  = 0.0;
            v4rho3sigma[12 * g + 7]  = 0.0;
            v4rho3sigma[12 * g + 8]  = 0.0;
            v4rho3sigma[12 * g + 9]  = 0.0;
            v4rho3sigma[12 * g + 10] = 0.0;
            v4rho3sigma[12 * g + 11] = 0.0;

            v4rho2sigma2[18 * g + 1]  = 0.0;
            v4rho2sigma2[18 * g + 2]  = 0.0;
            v4rho2sigma2[18 * g + 3]  = 0.0;
            v4rho2sigma2[18 * g + 4]  = 0.0;
            v4rho2sigma2[18 * g + 5]  = 0.0;
            v4rho2sigma2[18 * g + 6]  = 0.0;
            v4rho2sigma2[18 * g + 7]  = 0.0;
            v4rho2sigma2[18 * g + 8]  = 0.0;
            v4rho2sigma2[18 * g + 9]  = 0.0;
            v4rho2sigma2[18 * g + 10] = 0.0;
            v4rho2sigma2[18 * g + 11] = 0.0;
            v4rho2sigma2[18 * g + 12] = 0.0;
            v4rho2sigma2[18 * g + 13] = 0.0;
            v4rho2sigma2[18 * g + 14] = 0.0;
            v4rho2sigma2[18 * g + 15] = 0.0;
            v4rho2sigma2[18 * g + 16] = 0.0;
            v4rho2sigma2[18 * g + 17] = 0.0;

            v4rhosigma3[20 * g + 1]  = 0.0;
            v4rhosigma3[20 * g + 2]  = 0.0;
            v4rhosigma3[20 * g + 3]  = 0.0;
            v4rhosigma3[20 * g + 4]  = 0.0;
            v4rhosigma3[20 * g + 5]  = 0.0;
            v4rhosigma3[20 * g + 6]  = 0.0;
            v4rhosigma3[20 * g + 7]  = 0.0;
            v4rhosigma3[20 * g + 8]  = 0.0;
            v4rhosigma3[20 * g + 9]  = 0.0;
            v4rhosigma3[20 * g + 10] = 0.0;
            v4rhosigma3[20 * g + 11] = 0.0;
            v4rhosigma3[20 * g + 12] = 0.0;
            v4rhosigma3[20 * g + 13] = 0.0;
            v4rhosigma3[20 * g + 14] = 0.0;
            v4rhosigma3[20 * g + 15] = 0.0;
            v4rhosigma3[20 * g + 16] = 0.0;
            v4rhosigma3[20 * g + 17] = 0.0;
            v4rhosigma3[20 * g + 18] = 0.0;
            v4rhosigma3[20 * g + 19] = 0.0;

            v4sigma4[15 * g + 1]  = 0.0;
            v4sigma4[15 * g + 2]  = 0.0;
            v4sigma4[15 * g + 3]  = 0.0;
            v4sigma4[15 * g + 4]  = 0.0;
            v4sigma4[15 * g + 5]  = 0.0;
            v4sigma4[15 * g + 6]  = 0.0;
            v4sigma4[15 * g + 7]  = 0.0;
            v4sigma4[15 * g + 8]  = 0.0;
            v4sigma4[15 * g + 9]  = 0.0;
            v4sigma4[15 * g + 10] = 0.0;
            v4sigma4[15 * g + 11] = 0.0;
            v4sigma4[15 * g + 12] = 0.0;
            v4sigma4[15 * g + 13] = 0.0;
            v4sigma4[15 * g + 14] = 0.0;
        }
    }
}

void
screenLxcForMGGA(const int32_t npoints,
                 const double* rho,
                 const double* sigma,
                 const double* lapl,
                 const double* tau,
                 double*       v4rho4,
                 double*       v4rho3sigma,
                 double*       v4rho3lapl,
                 double*       v4rho3tau,
                 double*       v4rho2sigma2,
                 double*       v4rho2sigmalapl,
                 double*       v4rho2sigmatau,
                 double*       v4rho2lapl2,
                 double*       v4rho2lapltau,
                 double*       v4rho2tau2,
                 double*       v4rhosigma3,
                 double*       v4rhosigma2lapl,
                 double*       v4rhosigma2tau,
                 double*       v4rhosigmalapl2,
                 double*       v4rhosigmalapltau,
                 double*       v4rhosigmatau2,
                 double*       v4rholapl3,
                 double*       v4rholapl2tau,
                 double*       v4rholapltau2,
                 double*       v4rhotau3,
                 double*       v4sigma4,
                 double*       v4sigma3lapl,
                 double*       v4sigma3tau,
                 double*       v4sigma2lapl2,
                 double*       v4sigma2lapltau,
                 double*       v4sigma2tau2,
                 double*       v4sigmalapl3,
                 double*       v4sigmalapl2tau,
                 double*       v4sigmalapltau2,
                 double*       v4sigmatau3,
                 double*       v4lapl4,
                 double*       v4lapl3tau,
                 double*       v4lapl2tau2,
                 double*       v4lapltau3,
                 double*       v4tau4)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    for (int32_t g = 0; g < npoints; g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            v4rho4[5 * g + 0] = 0.0;
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;

            v4rho3sigma[12 * g + 0]  = 0.0;
            v4rho3sigma[12 * g + 1]  = 0.0;
            v4rho3sigma[12 * g + 2]  = 0.0;
            v4rho3sigma[12 * g + 3]  = 0.0;
            v4rho3sigma[12 * g + 4]  = 0.0;
            v4rho3sigma[12 * g + 5]  = 0.0;
            v4rho3sigma[12 * g + 6]  = 0.0;
            v4rho3sigma[12 * g + 7]  = 0.0;
            v4rho3sigma[12 * g + 8]  = 0.0;
            v4rho3sigma[12 * g + 9]  = 0.0;
            v4rho3sigma[12 * g + 10] = 0.0;

            v4rho3lapl[8 * g + 0] = 0.0;
            v4rho3lapl[8 * g + 1] = 0.0;
            v4rho3lapl[8 * g + 2] = 0.0;
            v4rho3lapl[8 * g + 3] = 0.0;
            v4rho3lapl[8 * g + 4] = 0.0;
            v4rho3lapl[8 * g + 5] = 0.0;
            v4rho3lapl[8 * g + 6] = 0.0;

            v4rho3tau[8 * g + 0] = 0.0;
            v4rho3tau[8 * g + 1] = 0.0;
            v4rho3tau[8 * g + 2] = 0.0;
            v4rho3tau[8 * g + 3] = 0.0;
            v4rho3tau[8 * g + 4] = 0.0;
            v4rho3tau[8 * g + 5] = 0.0;
            v4rho3tau[8 * g + 6] = 0.0;

            v4rho2sigma2[18 * g + 0]  = 0.0;
            v4rho2sigma2[18 * g + 1]  = 0.0;
            v4rho2sigma2[18 * g + 2]  = 0.0;
            v4rho2sigma2[18 * g + 3]  = 0.0;
            v4rho2sigma2[18 * g + 4]  = 0.0;
            v4rho2sigma2[18 * g + 5]  = 0.0;
            v4rho2sigma2[18 * g + 6]  = 0.0;
            v4rho2sigma2[18 * g + 7]  = 0.0;
            v4rho2sigma2[18 * g + 8]  = 0.0;
            v4rho2sigma2[18 * g + 9]  = 0.0;
            v4rho2sigma2[18 * g + 10] = 0.0;
            v4rho2sigma2[18 * g + 11] = 0.0;
            v4rho2sigma2[18 * g + 12] = 0.0;
            v4rho2sigma2[18 * g + 13] = 0.0;
            v4rho2sigma2[18 * g + 14] = 0.0;
            v4rho2sigma2[18 * g + 15] = 0.0;
            v4rho2sigma2[18 * g + 16] = 0.0;

            v4rho2sigmalapl[18 * g + 0]  = 0.0;
            v4rho2sigmalapl[18 * g + 1]  = 0.0;
            v4rho2sigmalapl[18 * g + 2]  = 0.0;
            v4rho2sigmalapl[18 * g + 3]  = 0.0;
            v4rho2sigmalapl[18 * g + 4]  = 0.0;
            v4rho2sigmalapl[18 * g + 5]  = 0.0;
            v4rho2sigmalapl[18 * g + 6]  = 0.0;
            v4rho2sigmalapl[18 * g + 7]  = 0.0;
            v4rho2sigmalapl[18 * g + 8]  = 0.0;
            v4rho2sigmalapl[18 * g + 9]  = 0.0;
            v4rho2sigmalapl[18 * g + 10] = 0.0;
            v4rho2sigmalapl[18 * g + 11] = 0.0;
            v4rho2sigmalapl[18 * g + 12] = 0.0;
            v4rho2sigmalapl[18 * g + 13] = 0.0;
            v4rho2sigmalapl[18 * g + 14] = 0.0;
            v4rho2sigmalapl[18 * g + 15] = 0.0;
            v4rho2sigmalapl[18 * g + 16] = 0.0;

            v4rho2sigmatau[18 * g + 0]  = 0.0;
            v4rho2sigmatau[18 * g + 1]  = 0.0;
            v4rho2sigmatau[18 * g + 2]  = 0.0;
            v4rho2sigmatau[18 * g + 3]  = 0.0;
            v4rho2sigmatau[18 * g + 4]  = 0.0;
            v4rho2sigmatau[18 * g + 5]  = 0.0;
            v4rho2sigmatau[18 * g + 6]  = 0.0;
            v4rho2sigmatau[18 * g + 7]  = 0.0;
            v4rho2sigmatau[18 * g + 8]  = 0.0;
            v4rho2sigmatau[18 * g + 9]  = 0.0;
            v4rho2sigmatau[18 * g + 10] = 0.0;
            v4rho2sigmatau[18 * g + 11] = 0.0;
            v4rho2sigmatau[18 * g + 12] = 0.0;
            v4rho2sigmatau[18 * g + 13] = 0.0;
            v4rho2sigmatau[18 * g + 14] = 0.0;
            v4rho2sigmatau[18 * g + 15] = 0.0;
            v4rho2sigmatau[18 * g + 16] = 0.0;

            v4rho2lapl2[9 * g + 0] = 0.0;
            v4rho2lapl2[9 * g + 1] = 0.0;
            v4rho2lapl2[9 * g + 2] = 0.0;
            v4rho2lapl2[9 * g + 3] = 0.0;
            v4rho2lapl2[9 * g + 4] = 0.0;
            v4rho2lapl2[9 * g + 5] = 0.0;
            v4rho2lapl2[9 * g + 6] = 0.0;
            v4rho2lapl2[9 * g + 7] = 0.0;

            v4rho2lapltau[12 * g + 0]  = 0.0;
            v4rho2lapltau[12 * g + 1]  = 0.0;
            v4rho2lapltau[12 * g + 2]  = 0.0;
            v4rho2lapltau[12 * g + 3]  = 0.0;
            v4rho2lapltau[12 * g + 4]  = 0.0;
            v4rho2lapltau[12 * g + 5]  = 0.0;
            v4rho2lapltau[12 * g + 6]  = 0.0;
            v4rho2lapltau[12 * g + 7]  = 0.0;
            v4rho2lapltau[12 * g + 8]  = 0.0;
            v4rho2lapltau[12 * g + 9]  = 0.0;
            v4rho2lapltau[12 * g + 10] = 0.0;

            v4rho2tau2[9 * g + 0] = 0.0;
            v4rho2tau2[9 * g + 1] = 0.0;
            v4rho2tau2[9 * g + 2] = 0.0;
            v4rho2tau2[9 * g + 3] = 0.0;
            v4rho2tau2[9 * g + 4] = 0.0;
            v4rho2tau2[9 * g + 5] = 0.0;
            v4rho2tau2[9 * g + 6] = 0.0;
            v4rho2tau2[9 * g + 7] = 0.0;

            v4rhosigma3[20 * g + 0]  = 0.0;
            v4rhosigma3[20 * g + 1]  = 0.0;
            v4rhosigma3[20 * g + 2]  = 0.0;
            v4rhosigma3[20 * g + 3]  = 0.0;
            v4rhosigma3[20 * g + 4]  = 0.0;
            v4rhosigma3[20 * g + 5]  = 0.0;
            v4rhosigma3[20 * g + 6]  = 0.0;
            v4rhosigma3[20 * g + 7]  = 0.0;
            v4rhosigma3[20 * g + 8]  = 0.0;
            v4rhosigma3[20 * g + 9]  = 0.0;
            v4rhosigma3[20 * g + 10] = 0.0;
            v4rhosigma3[20 * g + 11] = 0.0;
            v4rhosigma3[20 * g + 12] = 0.0;
            v4rhosigma3[20 * g + 13] = 0.0;
            v4rhosigma3[20 * g + 14] = 0.0;
            v4rhosigma3[20 * g + 15] = 0.0;
            v4rhosigma3[20 * g + 16] = 0.0;
            v4rhosigma3[20 * g + 17] = 0.0;
            v4rhosigma3[20 * g + 18] = 0.0;

            // v4rhosigma2lapl: inconsistent size in libxc (36 vs 24)
            v4rhosigma2lapl[36 * g + 0]  = 0.0;
            v4rhosigma2lapl[36 * g + 1]  = 0.0;
            v4rhosigma2lapl[36 * g + 2]  = 0.0;
            v4rhosigma2lapl[36 * g + 3]  = 0.0;
            v4rhosigma2lapl[36 * g + 4]  = 0.0;
            v4rhosigma2lapl[36 * g + 5]  = 0.0;
            v4rhosigma2lapl[36 * g + 6]  = 0.0;
            v4rhosigma2lapl[36 * g + 7]  = 0.0;
            v4rhosigma2lapl[36 * g + 8]  = 0.0;
            v4rhosigma2lapl[36 * g + 9]  = 0.0;
            v4rhosigma2lapl[36 * g + 10] = 0.0;
            v4rhosigma2lapl[36 * g + 11] = 0.0;
            v4rhosigma2lapl[36 * g + 12] = 0.0;
            v4rhosigma2lapl[36 * g + 13] = 0.0;
            v4rhosigma2lapl[36 * g + 14] = 0.0;
            v4rhosigma2lapl[36 * g + 15] = 0.0;
            v4rhosigma2lapl[36 * g + 16] = 0.0;
            v4rhosigma2lapl[36 * g + 17] = 0.0;
            v4rhosigma2lapl[36 * g + 18] = 0.0;
            v4rhosigma2lapl[36 * g + 19] = 0.0;
            v4rhosigma2lapl[36 * g + 20] = 0.0;
            v4rhosigma2lapl[36 * g + 21] = 0.0;
            v4rhosigma2lapl[36 * g + 22] = 0.0;

            // v4rhosigma2tau: inconsistent size in libxc (36 vs 24)
            v4rhosigma2tau[36 * g + 0]  = 0.0;
            v4rhosigma2tau[36 * g + 1]  = 0.0;
            v4rhosigma2tau[36 * g + 2]  = 0.0;
            v4rhosigma2tau[36 * g + 3]  = 0.0;
            v4rhosigma2tau[36 * g + 4]  = 0.0;
            v4rhosigma2tau[36 * g + 5]  = 0.0;
            v4rhosigma2tau[36 * g + 6]  = 0.0;
            v4rhosigma2tau[36 * g + 7]  = 0.0;
            v4rhosigma2tau[36 * g + 8]  = 0.0;
            v4rhosigma2tau[36 * g + 9]  = 0.0;
            v4rhosigma2tau[36 * g + 10] = 0.0;
            v4rhosigma2tau[36 * g + 11] = 0.0;
            v4rhosigma2tau[36 * g + 12] = 0.0;
            v4rhosigma2tau[36 * g + 13] = 0.0;
            v4rhosigma2tau[36 * g + 14] = 0.0;
            v4rhosigma2tau[36 * g + 15] = 0.0;
            v4rhosigma2tau[36 * g + 16] = 0.0;
            v4rhosigma2tau[36 * g + 17] = 0.0;
            v4rhosigma2tau[36 * g + 18] = 0.0;
            v4rhosigma2tau[36 * g + 19] = 0.0;
            v4rhosigma2tau[36 * g + 20] = 0.0;
            v4rhosigma2tau[36 * g + 21] = 0.0;
            v4rhosigma2tau[36 * g + 22] = 0.0;

            v4rhosigmalapl2[18 * g + 0]  = 0.0;
            v4rhosigmalapl2[18 * g + 1]  = 0.0;
            v4rhosigmalapl2[18 * g + 2]  = 0.0;
            v4rhosigmalapl2[18 * g + 3]  = 0.0;
            v4rhosigmalapl2[18 * g + 4]  = 0.0;
            v4rhosigmalapl2[18 * g + 5]  = 0.0;
            v4rhosigmalapl2[18 * g + 6]  = 0.0;
            v4rhosigmalapl2[18 * g + 7]  = 0.0;
            v4rhosigmalapl2[18 * g + 8]  = 0.0;
            v4rhosigmalapl2[18 * g + 9]  = 0.0;
            v4rhosigmalapl2[18 * g + 10] = 0.0;
            v4rhosigmalapl2[18 * g + 11] = 0.0;
            v4rhosigmalapl2[18 * g + 12] = 0.0;
            v4rhosigmalapl2[18 * g + 13] = 0.0;
            v4rhosigmalapl2[18 * g + 14] = 0.0;
            v4rhosigmalapl2[18 * g + 15] = 0.0;
            v4rhosigmalapl2[18 * g + 16] = 0.0;

            v4rhosigmalapltau[24 * g + 0]  = 0.0;
            v4rhosigmalapltau[24 * g + 1]  = 0.0;
            v4rhosigmalapltau[24 * g + 2]  = 0.0;
            v4rhosigmalapltau[24 * g + 3]  = 0.0;
            v4rhosigmalapltau[24 * g + 4]  = 0.0;
            v4rhosigmalapltau[24 * g + 5]  = 0.0;
            v4rhosigmalapltau[24 * g + 6]  = 0.0;
            v4rhosigmalapltau[24 * g + 7]  = 0.0;
            v4rhosigmalapltau[24 * g + 8]  = 0.0;
            v4rhosigmalapltau[24 * g + 9]  = 0.0;
            v4rhosigmalapltau[24 * g + 10] = 0.0;
            v4rhosigmalapltau[24 * g + 11] = 0.0;
            v4rhosigmalapltau[24 * g + 12] = 0.0;
            v4rhosigmalapltau[24 * g + 13] = 0.0;
            v4rhosigmalapltau[24 * g + 14] = 0.0;
            v4rhosigmalapltau[24 * g + 15] = 0.0;
            v4rhosigmalapltau[24 * g + 16] = 0.0;
            v4rhosigmalapltau[24 * g + 17] = 0.0;
            v4rhosigmalapltau[24 * g + 18] = 0.0;
            v4rhosigmalapltau[24 * g + 19] = 0.0;
            v4rhosigmalapltau[24 * g + 20] = 0.0;
            v4rhosigmalapltau[24 * g + 21] = 0.0;
            v4rhosigmalapltau[24 * g + 22] = 0.0;

            // v4rhosigmatau2: inconsistent size in libxc (36 vs 18)
            v4rhosigmatau2[36 * g + 0]  = 0.0;
            v4rhosigmatau2[36 * g + 1]  = 0.0;
            v4rhosigmatau2[36 * g + 2]  = 0.0;
            v4rhosigmatau2[36 * g + 3]  = 0.0;
            v4rhosigmatau2[36 * g + 4]  = 0.0;
            v4rhosigmatau2[36 * g + 5]  = 0.0;
            v4rhosigmatau2[36 * g + 6]  = 0.0;
            v4rhosigmatau2[36 * g + 7]  = 0.0;
            v4rhosigmatau2[36 * g + 8]  = 0.0;
            v4rhosigmatau2[36 * g + 9]  = 0.0;
            v4rhosigmatau2[36 * g + 10] = 0.0;
            v4rhosigmatau2[36 * g + 11] = 0.0;
            v4rhosigmatau2[36 * g + 12] = 0.0;
            v4rhosigmatau2[36 * g + 13] = 0.0;
            v4rhosigmatau2[36 * g + 14] = 0.0;
            v4rhosigmatau2[36 * g + 15] = 0.0;
            v4rhosigmatau2[36 * g + 16] = 0.0;

            v4rholapl3[8 * g + 0] = 0.0;
            v4rholapl3[8 * g + 1] = 0.0;
            v4rholapl3[8 * g + 2] = 0.0;
            v4rholapl3[8 * g + 3] = 0.0;
            v4rholapl3[8 * g + 4] = 0.0;
            v4rholapl3[8 * g + 5] = 0.0;
            v4rholapl3[8 * g + 6] = 0.0;

            v4rholapl2tau[12 * g + 0]  = 0.0;
            v4rholapl2tau[12 * g + 1]  = 0.0;
            v4rholapl2tau[12 * g + 2]  = 0.0;
            v4rholapl2tau[12 * g + 3]  = 0.0;
            v4rholapl2tau[12 * g + 4]  = 0.0;
            v4rholapl2tau[12 * g + 5]  = 0.0;
            v4rholapl2tau[12 * g + 6]  = 0.0;
            v4rholapl2tau[12 * g + 7]  = 0.0;
            v4rholapl2tau[12 * g + 8]  = 0.0;
            v4rholapl2tau[12 * g + 9]  = 0.0;
            v4rholapl2tau[12 * g + 10] = 0.0;

            v4rholapltau2[12 * g + 0]  = 0.0;
            v4rholapltau2[12 * g + 1]  = 0.0;
            v4rholapltau2[12 * g + 2]  = 0.0;
            v4rholapltau2[12 * g + 3]  = 0.0;
            v4rholapltau2[12 * g + 4]  = 0.0;
            v4rholapltau2[12 * g + 5]  = 0.0;
            v4rholapltau2[12 * g + 6]  = 0.0;
            v4rholapltau2[12 * g + 7]  = 0.0;
            v4rholapltau2[12 * g + 8]  = 0.0;
            v4rholapltau2[12 * g + 9]  = 0.0;
            v4rholapltau2[12 * g + 10] = 0.0;

            v4rhotau3[8 * g + 0] = 0.0;
            v4rhotau3[8 * g + 1] = 0.0;
            v4rhotau3[8 * g + 2] = 0.0;
            v4rhotau3[8 * g + 3] = 0.0;
            v4rhotau3[8 * g + 4] = 0.0;
            v4rhotau3[8 * g + 5] = 0.0;
            v4rhotau3[8 * g + 6] = 0.0;

            v4sigma4[15 * g + 0]  = 0.0;
            v4sigma4[15 * g + 1]  = 0.0;
            v4sigma4[15 * g + 2]  = 0.0;
            v4sigma4[15 * g + 3]  = 0.0;
            v4sigma4[15 * g + 4]  = 0.0;
            v4sigma4[15 * g + 5]  = 0.0;
            v4sigma4[15 * g + 6]  = 0.0;
            v4sigma4[15 * g + 7]  = 0.0;
            v4sigma4[15 * g + 8]  = 0.0;
            v4sigma4[15 * g + 9]  = 0.0;
            v4sigma4[15 * g + 10] = 0.0;
            v4sigma4[15 * g + 11] = 0.0;
            v4sigma4[15 * g + 12] = 0.0;
            v4sigma4[15 * g + 13] = 0.0;

            v4sigma3lapl[20 * g + 0]  = 0.0;
            v4sigma3lapl[20 * g + 1]  = 0.0;
            v4sigma3lapl[20 * g + 2]  = 0.0;
            v4sigma3lapl[20 * g + 3]  = 0.0;
            v4sigma3lapl[20 * g + 4]  = 0.0;
            v4sigma3lapl[20 * g + 5]  = 0.0;
            v4sigma3lapl[20 * g + 6]  = 0.0;
            v4sigma3lapl[20 * g + 7]  = 0.0;
            v4sigma3lapl[20 * g + 8]  = 0.0;
            v4sigma3lapl[20 * g + 9]  = 0.0;
            v4sigma3lapl[20 * g + 10] = 0.0;
            v4sigma3lapl[20 * g + 11] = 0.0;
            v4sigma3lapl[20 * g + 12] = 0.0;
            v4sigma3lapl[20 * g + 13] = 0.0;
            v4sigma3lapl[20 * g + 14] = 0.0;
            v4sigma3lapl[20 * g + 15] = 0.0;
            v4sigma3lapl[20 * g + 16] = 0.0;
            v4sigma3lapl[20 * g + 17] = 0.0;
            v4sigma3lapl[20 * g + 18] = 0.0;

            // v4sigma3tau: inconsistent size in libxc (30 vs 20)
            v4sigma3tau[30 * g + 0]  = 0.0;
            v4sigma3tau[30 * g + 1]  = 0.0;
            v4sigma3tau[30 * g + 2]  = 0.0;
            v4sigma3tau[30 * g + 3]  = 0.0;
            v4sigma3tau[30 * g + 4]  = 0.0;
            v4sigma3tau[30 * g + 5]  = 0.0;
            v4sigma3tau[30 * g + 6]  = 0.0;
            v4sigma3tau[30 * g + 7]  = 0.0;
            v4sigma3tau[30 * g + 8]  = 0.0;
            v4sigma3tau[30 * g + 9]  = 0.0;
            v4sigma3tau[30 * g + 10] = 0.0;
            v4sigma3tau[30 * g + 11] = 0.0;
            v4sigma3tau[30 * g + 12] = 0.0;
            v4sigma3tau[30 * g + 13] = 0.0;
            v4sigma3tau[30 * g + 14] = 0.0;
            v4sigma3tau[30 * g + 15] = 0.0;
            v4sigma3tau[30 * g + 16] = 0.0;
            v4sigma3tau[30 * g + 17] = 0.0;
            v4sigma3tau[30 * g + 18] = 0.0;

            v4sigma2lapl2[18 * g + 0]  = 0.0;
            v4sigma2lapl2[18 * g + 1]  = 0.0;
            v4sigma2lapl2[18 * g + 2]  = 0.0;
            v4sigma2lapl2[18 * g + 3]  = 0.0;
            v4sigma2lapl2[18 * g + 4]  = 0.0;
            v4sigma2lapl2[18 * g + 5]  = 0.0;
            v4sigma2lapl2[18 * g + 6]  = 0.0;
            v4sigma2lapl2[18 * g + 7]  = 0.0;
            v4sigma2lapl2[18 * g + 8]  = 0.0;
            v4sigma2lapl2[18 * g + 9]  = 0.0;
            v4sigma2lapl2[18 * g + 10] = 0.0;
            v4sigma2lapl2[18 * g + 11] = 0.0;
            v4sigma2lapl2[18 * g + 12] = 0.0;
            v4sigma2lapl2[18 * g + 13] = 0.0;
            v4sigma2lapl2[18 * g + 14] = 0.0;
            v4sigma2lapl2[18 * g + 15] = 0.0;
            v4sigma2lapl2[18 * g + 16] = 0.0;

            v4sigma2lapltau[24 * g + 0]  = 0.0;
            v4sigma2lapltau[24 * g + 1]  = 0.0;
            v4sigma2lapltau[24 * g + 2]  = 0.0;
            v4sigma2lapltau[24 * g + 3]  = 0.0;
            v4sigma2lapltau[24 * g + 4]  = 0.0;
            v4sigma2lapltau[24 * g + 5]  = 0.0;
            v4sigma2lapltau[24 * g + 6]  = 0.0;
            v4sigma2lapltau[24 * g + 7]  = 0.0;
            v4sigma2lapltau[24 * g + 8]  = 0.0;
            v4sigma2lapltau[24 * g + 9]  = 0.0;
            v4sigma2lapltau[24 * g + 10] = 0.0;
            v4sigma2lapltau[24 * g + 11] = 0.0;
            v4sigma2lapltau[24 * g + 12] = 0.0;
            v4sigma2lapltau[24 * g + 13] = 0.0;
            v4sigma2lapltau[24 * g + 14] = 0.0;
            v4sigma2lapltau[24 * g + 15] = 0.0;
            v4sigma2lapltau[24 * g + 16] = 0.0;
            v4sigma2lapltau[24 * g + 17] = 0.0;
            v4sigma2lapltau[24 * g + 18] = 0.0;
            v4sigma2lapltau[24 * g + 19] = 0.0;
            v4sigma2lapltau[24 * g + 20] = 0.0;
            v4sigma2lapltau[24 * g + 21] = 0.0;
            v4sigma2lapltau[24 * g + 22] = 0.0;

            v4sigma2tau2[18 * g + 0]  = 0.0;
            v4sigma2tau2[18 * g + 1]  = 0.0;
            v4sigma2tau2[18 * g + 2]  = 0.0;
            v4sigma2tau2[18 * g + 3]  = 0.0;
            v4sigma2tau2[18 * g + 4]  = 0.0;
            v4sigma2tau2[18 * g + 5]  = 0.0;
            v4sigma2tau2[18 * g + 6]  = 0.0;
            v4sigma2tau2[18 * g + 7]  = 0.0;
            v4sigma2tau2[18 * g + 8]  = 0.0;
            v4sigma2tau2[18 * g + 9]  = 0.0;
            v4sigma2tau2[18 * g + 10] = 0.0;
            v4sigma2tau2[18 * g + 11] = 0.0;
            v4sigma2tau2[18 * g + 12] = 0.0;
            v4sigma2tau2[18 * g + 13] = 0.0;
            v4sigma2tau2[18 * g + 14] = 0.0;
            v4sigma2tau2[18 * g + 15] = 0.0;
            v4sigma2tau2[18 * g + 16] = 0.0;

            v4sigmalapl3[12 * g + 0]  = 0.0;
            v4sigmalapl3[12 * g + 1]  = 0.0;
            v4sigmalapl3[12 * g + 2]  = 0.0;
            v4sigmalapl3[12 * g + 3]  = 0.0;
            v4sigmalapl3[12 * g + 4]  = 0.0;
            v4sigmalapl3[12 * g + 5]  = 0.0;
            v4sigmalapl3[12 * g + 6]  = 0.0;
            v4sigmalapl3[12 * g + 7]  = 0.0;
            v4sigmalapl3[12 * g + 8]  = 0.0;
            v4sigmalapl3[12 * g + 9]  = 0.0;
            v4sigmalapl3[12 * g + 10] = 0.0;

            v4sigmalapl2tau[18 * g + 0]  = 0.0;
            v4sigmalapl2tau[18 * g + 1]  = 0.0;
            v4sigmalapl2tau[18 * g + 2]  = 0.0;
            v4sigmalapl2tau[18 * g + 3]  = 0.0;
            v4sigmalapl2tau[18 * g + 4]  = 0.0;
            v4sigmalapl2tau[18 * g + 5]  = 0.0;
            v4sigmalapl2tau[18 * g + 6]  = 0.0;
            v4sigmalapl2tau[18 * g + 7]  = 0.0;
            v4sigmalapl2tau[18 * g + 8]  = 0.0;
            v4sigmalapl2tau[18 * g + 9]  = 0.0;
            v4sigmalapl2tau[18 * g + 10] = 0.0;
            v4sigmalapl2tau[18 * g + 11] = 0.0;
            v4sigmalapl2tau[18 * g + 12] = 0.0;
            v4sigmalapl2tau[18 * g + 13] = 0.0;
            v4sigmalapl2tau[18 * g + 14] = 0.0;
            v4sigmalapl2tau[18 * g + 15] = 0.0;
            v4sigmalapl2tau[18 * g + 16] = 0.0;

            v4sigmalapltau2[18 * g + 0]  = 0.0;
            v4sigmalapltau2[18 * g + 1]  = 0.0;
            v4sigmalapltau2[18 * g + 2]  = 0.0;
            v4sigmalapltau2[18 * g + 3]  = 0.0;
            v4sigmalapltau2[18 * g + 4]  = 0.0;
            v4sigmalapltau2[18 * g + 5]  = 0.0;
            v4sigmalapltau2[18 * g + 6]  = 0.0;
            v4sigmalapltau2[18 * g + 7]  = 0.0;
            v4sigmalapltau2[18 * g + 8]  = 0.0;
            v4sigmalapltau2[18 * g + 9]  = 0.0;
            v4sigmalapltau2[18 * g + 10] = 0.0;
            v4sigmalapltau2[18 * g + 11] = 0.0;
            v4sigmalapltau2[18 * g + 12] = 0.0;
            v4sigmalapltau2[18 * g + 13] = 0.0;
            v4sigmalapltau2[18 * g + 14] = 0.0;
            v4sigmalapltau2[18 * g + 15] = 0.0;
            v4sigmalapltau2[18 * g + 16] = 0.0;

            v4sigmatau3[12 * g + 0]  = 0.0;
            v4sigmatau3[12 * g + 1]  = 0.0;
            v4sigmatau3[12 * g + 2]  = 0.0;
            v4sigmatau3[12 * g + 3]  = 0.0;
            v4sigmatau3[12 * g + 4]  = 0.0;
            v4sigmatau3[12 * g + 5]  = 0.0;
            v4sigmatau3[12 * g + 6]  = 0.0;
            v4sigmatau3[12 * g + 7]  = 0.0;
            v4sigmatau3[12 * g + 8]  = 0.0;
            v4sigmatau3[12 * g + 9]  = 0.0;
            v4sigmatau3[12 * g + 10] = 0.0;

            v4lapl4[5 * g + 0] = 0.0;
            v4lapl4[5 * g + 1] = 0.0;
            v4lapl4[5 * g + 2] = 0.0;
            v4lapl4[5 * g + 3] = 0.0;

            v4lapl3tau[8 * g + 0] = 0.0;
            v4lapl3tau[8 * g + 1] = 0.0;
            v4lapl3tau[8 * g + 2] = 0.0;
            v4lapl3tau[8 * g + 3] = 0.0;
            v4lapl3tau[8 * g + 4] = 0.0;
            v4lapl3tau[8 * g + 5] = 0.0;
            v4lapl3tau[8 * g + 6] = 0.0;

            v4lapl2tau2[9 * g + 0] = 0.0;
            v4lapl2tau2[9 * g + 1] = 0.0;
            v4lapl2tau2[9 * g + 2] = 0.0;
            v4lapl2tau2[9 * g + 3] = 0.0;
            v4lapl2tau2[9 * g + 4] = 0.0;
            v4lapl2tau2[9 * g + 5] = 0.0;
            v4lapl2tau2[9 * g + 6] = 0.0;
            v4lapl2tau2[9 * g + 7] = 0.0;

            v4lapltau3[8 * g + 0] = 0.0;
            v4lapltau3[8 * g + 1] = 0.0;
            v4lapltau3[8 * g + 2] = 0.0;
            v4lapltau3[8 * g + 3] = 0.0;
            v4lapltau3[8 * g + 4] = 0.0;
            v4lapltau3[8 * g + 5] = 0.0;
            v4lapltau3[8 * g + 6] = 0.0;

            v4tau4[5 * g + 0] = 0.0;
            v4tau4[5 * g + 1] = 0.0;
            v4tau4[5 * g + 2] = 0.0;
            v4tau4[5 * g + 3] = 0.0;
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            v4rho4[5 * g + 1] = 0.0;
            v4rho4[5 * g + 2] = 0.0;
            v4rho4[5 * g + 3] = 0.0;
            v4rho4[5 * g + 4] = 0.0;

            v4rho3sigma[12 * g + 1]  = 0.0;
            v4rho3sigma[12 * g + 2]  = 0.0;
            v4rho3sigma[12 * g + 3]  = 0.0;
            v4rho3sigma[12 * g + 4]  = 0.0;
            v4rho3sigma[12 * g + 5]  = 0.0;
            v4rho3sigma[12 * g + 6]  = 0.0;
            v4rho3sigma[12 * g + 7]  = 0.0;
            v4rho3sigma[12 * g + 8]  = 0.0;
            v4rho3sigma[12 * g + 9]  = 0.0;
            v4rho3sigma[12 * g + 10] = 0.0;
            v4rho3sigma[12 * g + 11] = 0.0;

            v4rho3lapl[8 * g + 1] = 0.0;
            v4rho3lapl[8 * g + 2] = 0.0;
            v4rho3lapl[8 * g + 3] = 0.0;
            v4rho3lapl[8 * g + 4] = 0.0;
            v4rho3lapl[8 * g + 5] = 0.0;
            v4rho3lapl[8 * g + 6] = 0.0;
            v4rho3lapl[8 * g + 7] = 0.0;

            v4rho3tau[8 * g + 1] = 0.0;
            v4rho3tau[8 * g + 2] = 0.0;
            v4rho3tau[8 * g + 3] = 0.0;
            v4rho3tau[8 * g + 4] = 0.0;
            v4rho3tau[8 * g + 5] = 0.0;
            v4rho3tau[8 * g + 6] = 0.0;
            v4rho3tau[8 * g + 7] = 0.0;

            v4rho2sigma2[18 * g + 1]  = 0.0;
            v4rho2sigma2[18 * g + 2]  = 0.0;
            v4rho2sigma2[18 * g + 3]  = 0.0;
            v4rho2sigma2[18 * g + 4]  = 0.0;
            v4rho2sigma2[18 * g + 5]  = 0.0;
            v4rho2sigma2[18 * g + 6]  = 0.0;
            v4rho2sigma2[18 * g + 7]  = 0.0;
            v4rho2sigma2[18 * g + 8]  = 0.0;
            v4rho2sigma2[18 * g + 9]  = 0.0;
            v4rho2sigma2[18 * g + 10] = 0.0;
            v4rho2sigma2[18 * g + 11] = 0.0;
            v4rho2sigma2[18 * g + 12] = 0.0;
            v4rho2sigma2[18 * g + 13] = 0.0;
            v4rho2sigma2[18 * g + 14] = 0.0;
            v4rho2sigma2[18 * g + 15] = 0.0;
            v4rho2sigma2[18 * g + 16] = 0.0;
            v4rho2sigma2[18 * g + 17] = 0.0;

            v4rho2sigmalapl[18 * g + 1]  = 0.0;
            v4rho2sigmalapl[18 * g + 2]  = 0.0;
            v4rho2sigmalapl[18 * g + 3]  = 0.0;
            v4rho2sigmalapl[18 * g + 4]  = 0.0;
            v4rho2sigmalapl[18 * g + 5]  = 0.0;
            v4rho2sigmalapl[18 * g + 6]  = 0.0;
            v4rho2sigmalapl[18 * g + 7]  = 0.0;
            v4rho2sigmalapl[18 * g + 8]  = 0.0;
            v4rho2sigmalapl[18 * g + 9]  = 0.0;
            v4rho2sigmalapl[18 * g + 10] = 0.0;
            v4rho2sigmalapl[18 * g + 11] = 0.0;
            v4rho2sigmalapl[18 * g + 12] = 0.0;
            v4rho2sigmalapl[18 * g + 13] = 0.0;
            v4rho2sigmalapl[18 * g + 14] = 0.0;
            v4rho2sigmalapl[18 * g + 15] = 0.0;
            v4rho2sigmalapl[18 * g + 16] = 0.0;
            v4rho2sigmalapl[18 * g + 17] = 0.0;

            v4rho2sigmatau[18 * g + 1]  = 0.0;
            v4rho2sigmatau[18 * g + 2]  = 0.0;
            v4rho2sigmatau[18 * g + 3]  = 0.0;
            v4rho2sigmatau[18 * g + 4]  = 0.0;
            v4rho2sigmatau[18 * g + 5]  = 0.0;
            v4rho2sigmatau[18 * g + 6]  = 0.0;
            v4rho2sigmatau[18 * g + 7]  = 0.0;
            v4rho2sigmatau[18 * g + 8]  = 0.0;
            v4rho2sigmatau[18 * g + 9]  = 0.0;
            v4rho2sigmatau[18 * g + 10] = 0.0;
            v4rho2sigmatau[18 * g + 11] = 0.0;
            v4rho2sigmatau[18 * g + 12] = 0.0;
            v4rho2sigmatau[18 * g + 13] = 0.0;
            v4rho2sigmatau[18 * g + 14] = 0.0;
            v4rho2sigmatau[18 * g + 15] = 0.0;
            v4rho2sigmatau[18 * g + 16] = 0.0;
            v4rho2sigmatau[18 * g + 17] = 0.0;

            v4rho2lapl2[9 * g + 1] = 0.0;
            v4rho2lapl2[9 * g + 2] = 0.0;
            v4rho2lapl2[9 * g + 3] = 0.0;
            v4rho2lapl2[9 * g + 4] = 0.0;
            v4rho2lapl2[9 * g + 5] = 0.0;
            v4rho2lapl2[9 * g + 6] = 0.0;
            v4rho2lapl2[9 * g + 7] = 0.0;
            v4rho2lapl2[9 * g + 8] = 0.0;

            v4rho2lapltau[12 * g + 1]  = 0.0;
            v4rho2lapltau[12 * g + 2]  = 0.0;
            v4rho2lapltau[12 * g + 3]  = 0.0;
            v4rho2lapltau[12 * g + 4]  = 0.0;
            v4rho2lapltau[12 * g + 5]  = 0.0;
            v4rho2lapltau[12 * g + 6]  = 0.0;
            v4rho2lapltau[12 * g + 7]  = 0.0;
            v4rho2lapltau[12 * g + 8]  = 0.0;
            v4rho2lapltau[12 * g + 9]  = 0.0;
            v4rho2lapltau[12 * g + 10] = 0.0;
            v4rho2lapltau[12 * g + 11] = 0.0;

            v4rho2tau2[9 * g + 1] = 0.0;
            v4rho2tau2[9 * g + 2] = 0.0;
            v4rho2tau2[9 * g + 3] = 0.0;
            v4rho2tau2[9 * g + 4] = 0.0;
            v4rho2tau2[9 * g + 5] = 0.0;
            v4rho2tau2[9 * g + 6] = 0.0;
            v4rho2tau2[9 * g + 7] = 0.0;
            v4rho2tau2[9 * g + 8] = 0.0;

            v4rhosigma3[20 * g + 1]  = 0.0;
            v4rhosigma3[20 * g + 2]  = 0.0;
            v4rhosigma3[20 * g + 3]  = 0.0;
            v4rhosigma3[20 * g + 4]  = 0.0;
            v4rhosigma3[20 * g + 5]  = 0.0;
            v4rhosigma3[20 * g + 6]  = 0.0;
            v4rhosigma3[20 * g + 7]  = 0.0;
            v4rhosigma3[20 * g + 8]  = 0.0;
            v4rhosigma3[20 * g + 9]  = 0.0;
            v4rhosigma3[20 * g + 10] = 0.0;
            v4rhosigma3[20 * g + 11] = 0.0;
            v4rhosigma3[20 * g + 12] = 0.0;
            v4rhosigma3[20 * g + 13] = 0.0;
            v4rhosigma3[20 * g + 14] = 0.0;
            v4rhosigma3[20 * g + 15] = 0.0;
            v4rhosigma3[20 * g + 16] = 0.0;
            v4rhosigma3[20 * g + 17] = 0.0;
            v4rhosigma3[20 * g + 18] = 0.0;
            v4rhosigma3[20 * g + 19] = 0.0;

            // v4rhosigma2lapl: inconsistent size in libxc (36 vs 24)
            v4rhosigma2lapl[36 * g + 1]  = 0.0;
            v4rhosigma2lapl[36 * g + 2]  = 0.0;
            v4rhosigma2lapl[36 * g + 3]  = 0.0;
            v4rhosigma2lapl[36 * g + 4]  = 0.0;
            v4rhosigma2lapl[36 * g + 5]  = 0.0;
            v4rhosigma2lapl[36 * g + 6]  = 0.0;
            v4rhosigma2lapl[36 * g + 7]  = 0.0;
            v4rhosigma2lapl[36 * g + 8]  = 0.0;
            v4rhosigma2lapl[36 * g + 9]  = 0.0;
            v4rhosigma2lapl[36 * g + 10] = 0.0;
            v4rhosigma2lapl[36 * g + 11] = 0.0;
            v4rhosigma2lapl[36 * g + 12] = 0.0;
            v4rhosigma2lapl[36 * g + 13] = 0.0;
            v4rhosigma2lapl[36 * g + 14] = 0.0;
            v4rhosigma2lapl[36 * g + 15] = 0.0;
            v4rhosigma2lapl[36 * g + 16] = 0.0;
            v4rhosigma2lapl[36 * g + 17] = 0.0;
            v4rhosigma2lapl[36 * g + 18] = 0.0;
            v4rhosigma2lapl[36 * g + 19] = 0.0;
            v4rhosigma2lapl[36 * g + 20] = 0.0;
            v4rhosigma2lapl[36 * g + 21] = 0.0;
            v4rhosigma2lapl[36 * g + 22] = 0.0;
            v4rhosigma2lapl[36 * g + 23] = 0.0;

            // v4rhosigma2tau: inconsistent size in libxc (36 vs 24)
            v4rhosigma2tau[36 * g + 1]  = 0.0;
            v4rhosigma2tau[36 * g + 2]  = 0.0;
            v4rhosigma2tau[36 * g + 3]  = 0.0;
            v4rhosigma2tau[36 * g + 4]  = 0.0;
            v4rhosigma2tau[36 * g + 5]  = 0.0;
            v4rhosigma2tau[36 * g + 6]  = 0.0;
            v4rhosigma2tau[36 * g + 7]  = 0.0;
            v4rhosigma2tau[36 * g + 8]  = 0.0;
            v4rhosigma2tau[36 * g + 9]  = 0.0;
            v4rhosigma2tau[36 * g + 10] = 0.0;
            v4rhosigma2tau[36 * g + 11] = 0.0;
            v4rhosigma2tau[36 * g + 12] = 0.0;
            v4rhosigma2tau[36 * g + 13] = 0.0;
            v4rhosigma2tau[36 * g + 14] = 0.0;
            v4rhosigma2tau[36 * g + 15] = 0.0;
            v4rhosigma2tau[36 * g + 16] = 0.0;
            v4rhosigma2tau[36 * g + 17] = 0.0;
            v4rhosigma2tau[36 * g + 18] = 0.0;
            v4rhosigma2tau[36 * g + 19] = 0.0;
            v4rhosigma2tau[36 * g + 20] = 0.0;
            v4rhosigma2tau[36 * g + 21] = 0.0;
            v4rhosigma2tau[36 * g + 22] = 0.0;
            v4rhosigma2tau[36 * g + 23] = 0.0;

            v4rhosigmalapl2[18 * g + 1]  = 0.0;
            v4rhosigmalapl2[18 * g + 2]  = 0.0;
            v4rhosigmalapl2[18 * g + 3]  = 0.0;
            v4rhosigmalapl2[18 * g + 4]  = 0.0;
            v4rhosigmalapl2[18 * g + 5]  = 0.0;
            v4rhosigmalapl2[18 * g + 6]  = 0.0;
            v4rhosigmalapl2[18 * g + 7]  = 0.0;
            v4rhosigmalapl2[18 * g + 8]  = 0.0;
            v4rhosigmalapl2[18 * g + 9]  = 0.0;
            v4rhosigmalapl2[18 * g + 10] = 0.0;
            v4rhosigmalapl2[18 * g + 11] = 0.0;
            v4rhosigmalapl2[18 * g + 12] = 0.0;
            v4rhosigmalapl2[18 * g + 13] = 0.0;
            v4rhosigmalapl2[18 * g + 14] = 0.0;
            v4rhosigmalapl2[18 * g + 15] = 0.0;
            v4rhosigmalapl2[18 * g + 16] = 0.0;
            v4rhosigmalapl2[18 * g + 17] = 0.0;

            v4rhosigmalapltau[24 * g + 1]  = 0.0;
            v4rhosigmalapltau[24 * g + 2]  = 0.0;
            v4rhosigmalapltau[24 * g + 3]  = 0.0;
            v4rhosigmalapltau[24 * g + 4]  = 0.0;
            v4rhosigmalapltau[24 * g + 5]  = 0.0;
            v4rhosigmalapltau[24 * g + 6]  = 0.0;
            v4rhosigmalapltau[24 * g + 7]  = 0.0;
            v4rhosigmalapltau[24 * g + 8]  = 0.0;
            v4rhosigmalapltau[24 * g + 9]  = 0.0;
            v4rhosigmalapltau[24 * g + 10] = 0.0;
            v4rhosigmalapltau[24 * g + 11] = 0.0;
            v4rhosigmalapltau[24 * g + 12] = 0.0;
            v4rhosigmalapltau[24 * g + 13] = 0.0;
            v4rhosigmalapltau[24 * g + 14] = 0.0;
            v4rhosigmalapltau[24 * g + 15] = 0.0;
            v4rhosigmalapltau[24 * g + 16] = 0.0;
            v4rhosigmalapltau[24 * g + 17] = 0.0;
            v4rhosigmalapltau[24 * g + 18] = 0.0;
            v4rhosigmalapltau[24 * g + 19] = 0.0;
            v4rhosigmalapltau[24 * g + 20] = 0.0;
            v4rhosigmalapltau[24 * g + 21] = 0.0;
            v4rhosigmalapltau[24 * g + 22] = 0.0;
            v4rhosigmalapltau[24 * g + 23] = 0.0;

            // v4rhosigmatau2: inconsistent size in libxc (36 vs 18)
            v4rhosigmatau2[36 * g + 1]  = 0.0;
            v4rhosigmatau2[36 * g + 2]  = 0.0;
            v4rhosigmatau2[36 * g + 3]  = 0.0;
            v4rhosigmatau2[36 * g + 4]  = 0.0;
            v4rhosigmatau2[36 * g + 5]  = 0.0;
            v4rhosigmatau2[36 * g + 6]  = 0.0;
            v4rhosigmatau2[36 * g + 7]  = 0.0;
            v4rhosigmatau2[36 * g + 8]  = 0.0;
            v4rhosigmatau2[36 * g + 9]  = 0.0;
            v4rhosigmatau2[36 * g + 10] = 0.0;
            v4rhosigmatau2[36 * g + 11] = 0.0;
            v4rhosigmatau2[36 * g + 12] = 0.0;
            v4rhosigmatau2[36 * g + 13] = 0.0;
            v4rhosigmatau2[36 * g + 14] = 0.0;
            v4rhosigmatau2[36 * g + 15] = 0.0;
            v4rhosigmatau2[36 * g + 16] = 0.0;
            v4rhosigmatau2[36 * g + 17] = 0.0;

            v4rholapl3[8 * g + 1] = 0.0;
            v4rholapl3[8 * g + 2] = 0.0;
            v4rholapl3[8 * g + 3] = 0.0;
            v4rholapl3[8 * g + 4] = 0.0;
            v4rholapl3[8 * g + 5] = 0.0;
            v4rholapl3[8 * g + 6] = 0.0;
            v4rholapl3[8 * g + 7] = 0.0;

            v4rholapl2tau[12 * g + 1]  = 0.0;
            v4rholapl2tau[12 * g + 2]  = 0.0;
            v4rholapl2tau[12 * g + 3]  = 0.0;
            v4rholapl2tau[12 * g + 4]  = 0.0;
            v4rholapl2tau[12 * g + 5]  = 0.0;
            v4rholapl2tau[12 * g + 6]  = 0.0;
            v4rholapl2tau[12 * g + 7]  = 0.0;
            v4rholapl2tau[12 * g + 8]  = 0.0;
            v4rholapl2tau[12 * g + 9]  = 0.0;
            v4rholapl2tau[12 * g + 10] = 0.0;
            v4rholapl2tau[12 * g + 11] = 0.0;

            v4rholapltau2[12 * g + 1]  = 0.0;
            v4rholapltau2[12 * g + 2]  = 0.0;
            v4rholapltau2[12 * g + 3]  = 0.0;
            v4rholapltau2[12 * g + 4]  = 0.0;
            v4rholapltau2[12 * g + 5]  = 0.0;
            v4rholapltau2[12 * g + 6]  = 0.0;
            v4rholapltau2[12 * g + 7]  = 0.0;
            v4rholapltau2[12 * g + 8]  = 0.0;
            v4rholapltau2[12 * g + 9]  = 0.0;
            v4rholapltau2[12 * g + 10] = 0.0;
            v4rholapltau2[12 * g + 11] = 0.0;

            v4rhotau3[8 * g + 1] = 0.0;
            v4rhotau3[8 * g + 2] = 0.0;
            v4rhotau3[8 * g + 3] = 0.0;
            v4rhotau3[8 * g + 4] = 0.0;
            v4rhotau3[8 * g + 5] = 0.0;
            v4rhotau3[8 * g + 6] = 0.0;
            v4rhotau3[8 * g + 7] = 0.0;

            v4sigma4[15 * g + 1]  = 0.0;
            v4sigma4[15 * g + 2]  = 0.0;
            v4sigma4[15 * g + 3]  = 0.0;
            v4sigma4[15 * g + 4]  = 0.0;
            v4sigma4[15 * g + 5]  = 0.0;
            v4sigma4[15 * g + 6]  = 0.0;
            v4sigma4[15 * g + 7]  = 0.0;
            v4sigma4[15 * g + 8]  = 0.0;
            v4sigma4[15 * g + 9]  = 0.0;
            v4sigma4[15 * g + 10] = 0.0;
            v4sigma4[15 * g + 11] = 0.0;
            v4sigma4[15 * g + 12] = 0.0;
            v4sigma4[15 * g + 13] = 0.0;
            v4sigma4[15 * g + 14] = 0.0;

            v4sigma3lapl[20 * g + 1]  = 0.0;
            v4sigma3lapl[20 * g + 2]  = 0.0;
            v4sigma3lapl[20 * g + 3]  = 0.0;
            v4sigma3lapl[20 * g + 4]  = 0.0;
            v4sigma3lapl[20 * g + 5]  = 0.0;
            v4sigma3lapl[20 * g + 6]  = 0.0;
            v4sigma3lapl[20 * g + 7]  = 0.0;
            v4sigma3lapl[20 * g + 8]  = 0.0;
            v4sigma3lapl[20 * g + 9]  = 0.0;
            v4sigma3lapl[20 * g + 10] = 0.0;
            v4sigma3lapl[20 * g + 11] = 0.0;
            v4sigma3lapl[20 * g + 12] = 0.0;
            v4sigma3lapl[20 * g + 13] = 0.0;
            v4sigma3lapl[20 * g + 14] = 0.0;
            v4sigma3lapl[20 * g + 15] = 0.0;
            v4sigma3lapl[20 * g + 16] = 0.0;
            v4sigma3lapl[20 * g + 17] = 0.0;
            v4sigma3lapl[20 * g + 18] = 0.0;
            v4sigma3lapl[20 * g + 19] = 0.0;

            // v4sigma3tau: inconsistent size in libxc (30 vs 20)
            v4sigma3tau[30 * g + 1]  = 0.0;
            v4sigma3tau[30 * g + 2]  = 0.0;
            v4sigma3tau[30 * g + 3]  = 0.0;
            v4sigma3tau[30 * g + 4]  = 0.0;
            v4sigma3tau[30 * g + 5]  = 0.0;
            v4sigma3tau[30 * g + 6]  = 0.0;
            v4sigma3tau[30 * g + 7]  = 0.0;
            v4sigma3tau[30 * g + 8]  = 0.0;
            v4sigma3tau[30 * g + 9]  = 0.0;
            v4sigma3tau[30 * g + 10] = 0.0;
            v4sigma3tau[30 * g + 11] = 0.0;
            v4sigma3tau[30 * g + 12] = 0.0;
            v4sigma3tau[30 * g + 13] = 0.0;
            v4sigma3tau[30 * g + 14] = 0.0;
            v4sigma3tau[30 * g + 15] = 0.0;
            v4sigma3tau[30 * g + 16] = 0.0;
            v4sigma3tau[30 * g + 17] = 0.0;
            v4sigma3tau[30 * g + 18] = 0.0;
            v4sigma3tau[30 * g + 19] = 0.0;

            v4sigma2lapl2[18 * g + 1]  = 0.0;
            v4sigma2lapl2[18 * g + 2]  = 0.0;
            v4sigma2lapl2[18 * g + 3]  = 0.0;
            v4sigma2lapl2[18 * g + 4]  = 0.0;
            v4sigma2lapl2[18 * g + 5]  = 0.0;
            v4sigma2lapl2[18 * g + 6]  = 0.0;
            v4sigma2lapl2[18 * g + 7]  = 0.0;
            v4sigma2lapl2[18 * g + 8]  = 0.0;
            v4sigma2lapl2[18 * g + 9]  = 0.0;
            v4sigma2lapl2[18 * g + 10] = 0.0;
            v4sigma2lapl2[18 * g + 11] = 0.0;
            v4sigma2lapl2[18 * g + 12] = 0.0;
            v4sigma2lapl2[18 * g + 13] = 0.0;
            v4sigma2lapl2[18 * g + 14] = 0.0;
            v4sigma2lapl2[18 * g + 15] = 0.0;
            v4sigma2lapl2[18 * g + 16] = 0.0;
            v4sigma2lapl2[18 * g + 17] = 0.0;

            v4sigma2lapltau[24 * g + 1]  = 0.0;
            v4sigma2lapltau[24 * g + 2]  = 0.0;
            v4sigma2lapltau[24 * g + 3]  = 0.0;
            v4sigma2lapltau[24 * g + 4]  = 0.0;
            v4sigma2lapltau[24 * g + 5]  = 0.0;
            v4sigma2lapltau[24 * g + 6]  = 0.0;
            v4sigma2lapltau[24 * g + 7]  = 0.0;
            v4sigma2lapltau[24 * g + 8]  = 0.0;
            v4sigma2lapltau[24 * g + 9]  = 0.0;
            v4sigma2lapltau[24 * g + 10] = 0.0;
            v4sigma2lapltau[24 * g + 11] = 0.0;
            v4sigma2lapltau[24 * g + 12] = 0.0;
            v4sigma2lapltau[24 * g + 13] = 0.0;
            v4sigma2lapltau[24 * g + 14] = 0.0;
            v4sigma2lapltau[24 * g + 15] = 0.0;
            v4sigma2lapltau[24 * g + 16] = 0.0;
            v4sigma2lapltau[24 * g + 17] = 0.0;
            v4sigma2lapltau[24 * g + 18] = 0.0;
            v4sigma2lapltau[24 * g + 19] = 0.0;
            v4sigma2lapltau[24 * g + 20] = 0.0;
            v4sigma2lapltau[24 * g + 21] = 0.0;
            v4sigma2lapltau[24 * g + 22] = 0.0;
            v4sigma2lapltau[24 * g + 23] = 0.0;

            v4sigma2tau2[18 * g + 1]  = 0.0;
            v4sigma2tau2[18 * g + 2]  = 0.0;
            v4sigma2tau2[18 * g + 3]  = 0.0;
            v4sigma2tau2[18 * g + 4]  = 0.0;
            v4sigma2tau2[18 * g + 5]  = 0.0;
            v4sigma2tau2[18 * g + 6]  = 0.0;
            v4sigma2tau2[18 * g + 7]  = 0.0;
            v4sigma2tau2[18 * g + 8]  = 0.0;
            v4sigma2tau2[18 * g + 9]  = 0.0;
            v4sigma2tau2[18 * g + 10] = 0.0;
            v4sigma2tau2[18 * g + 11] = 0.0;
            v4sigma2tau2[18 * g + 12] = 0.0;
            v4sigma2tau2[18 * g + 13] = 0.0;
            v4sigma2tau2[18 * g + 14] = 0.0;
            v4sigma2tau2[18 * g + 15] = 0.0;
            v4sigma2tau2[18 * g + 16] = 0.0;
            v4sigma2tau2[18 * g + 17] = 0.0;

            v4sigmalapl3[12 * g + 1]  = 0.0;
            v4sigmalapl3[12 * g + 2]  = 0.0;
            v4sigmalapl3[12 * g + 3]  = 0.0;
            v4sigmalapl3[12 * g + 4]  = 0.0;
            v4sigmalapl3[12 * g + 5]  = 0.0;
            v4sigmalapl3[12 * g + 6]  = 0.0;
            v4sigmalapl3[12 * g + 7]  = 0.0;
            v4sigmalapl3[12 * g + 8]  = 0.0;
            v4sigmalapl3[12 * g + 9]  = 0.0;
            v4sigmalapl3[12 * g + 10] = 0.0;
            v4sigmalapl3[12 * g + 11] = 0.0;

            v4sigmalapl2tau[18 * g + 1]  = 0.0;
            v4sigmalapl2tau[18 * g + 2]  = 0.0;
            v4sigmalapl2tau[18 * g + 3]  = 0.0;
            v4sigmalapl2tau[18 * g + 4]  = 0.0;
            v4sigmalapl2tau[18 * g + 5]  = 0.0;
            v4sigmalapl2tau[18 * g + 6]  = 0.0;
            v4sigmalapl2tau[18 * g + 7]  = 0.0;
            v4sigmalapl2tau[18 * g + 8]  = 0.0;
            v4sigmalapl2tau[18 * g + 9]  = 0.0;
            v4sigmalapl2tau[18 * g + 10] = 0.0;
            v4sigmalapl2tau[18 * g + 11] = 0.0;
            v4sigmalapl2tau[18 * g + 12] = 0.0;
            v4sigmalapl2tau[18 * g + 13] = 0.0;
            v4sigmalapl2tau[18 * g + 14] = 0.0;
            v4sigmalapl2tau[18 * g + 15] = 0.0;
            v4sigmalapl2tau[18 * g + 16] = 0.0;
            v4sigmalapl2tau[18 * g + 17] = 0.0;

            v4sigmalapltau2[18 * g + 1]  = 0.0;
            v4sigmalapltau2[18 * g + 2]  = 0.0;
            v4sigmalapltau2[18 * g + 3]  = 0.0;
            v4sigmalapltau2[18 * g + 4]  = 0.0;
            v4sigmalapltau2[18 * g + 5]  = 0.0;
            v4sigmalapltau2[18 * g + 6]  = 0.0;
            v4sigmalapltau2[18 * g + 7]  = 0.0;
            v4sigmalapltau2[18 * g + 8]  = 0.0;
            v4sigmalapltau2[18 * g + 9]  = 0.0;
            v4sigmalapltau2[18 * g + 10] = 0.0;
            v4sigmalapltau2[18 * g + 11] = 0.0;
            v4sigmalapltau2[18 * g + 12] = 0.0;
            v4sigmalapltau2[18 * g + 13] = 0.0;
            v4sigmalapltau2[18 * g + 14] = 0.0;
            v4sigmalapltau2[18 * g + 15] = 0.0;
            v4sigmalapltau2[18 * g + 16] = 0.0;
            v4sigmalapltau2[18 * g + 17] = 0.0;

            v4sigmatau3[12 * g + 1]  = 0.0;
            v4sigmatau3[12 * g + 2]  = 0.0;
            v4sigmatau3[12 * g + 3]  = 0.0;
            v4sigmatau3[12 * g + 4]  = 0.0;
            v4sigmatau3[12 * g + 5]  = 0.0;
            v4sigmatau3[12 * g + 6]  = 0.0;
            v4sigmatau3[12 * g + 7]  = 0.0;
            v4sigmatau3[12 * g + 8]  = 0.0;
            v4sigmatau3[12 * g + 9]  = 0.0;
            v4sigmatau3[12 * g + 10] = 0.0;
            v4sigmatau3[12 * g + 11] = 0.0;

            v4lapl4[5 * g + 1] = 0.0;
            v4lapl4[5 * g + 2] = 0.0;
            v4lapl4[5 * g + 3] = 0.0;
            v4lapl4[5 * g + 4] = 0.0;

            v4lapl3tau[8 * g + 1] = 0.0;
            v4lapl3tau[8 * g + 2] = 0.0;
            v4lapl3tau[8 * g + 3] = 0.0;
            v4lapl3tau[8 * g + 4] = 0.0;
            v4lapl3tau[8 * g + 5] = 0.0;
            v4lapl3tau[8 * g + 6] = 0.0;
            v4lapl3tau[8 * g + 7] = 0.0;

            v4lapl2tau2[9 * g + 1] = 0.0;
            v4lapl2tau2[9 * g + 2] = 0.0;
            v4lapl2tau2[9 * g + 3] = 0.0;
            v4lapl2tau2[9 * g + 4] = 0.0;
            v4lapl2tau2[9 * g + 5] = 0.0;
            v4lapl2tau2[9 * g + 6] = 0.0;
            v4lapl2tau2[9 * g + 7] = 0.0;
            v4lapl2tau2[9 * g + 8] = 0.0;

            v4lapltau3[8 * g + 1] = 0.0;
            v4lapltau3[8 * g + 2] = 0.0;
            v4lapltau3[8 * g + 3] = 0.0;
            v4lapltau3[8 * g + 4] = 0.0;
            v4lapltau3[8 * g + 5] = 0.0;
            v4lapltau3[8 * g + 6] = 0.0;
            v4lapltau3[8 * g + 7] = 0.0;

            v4tau4[5 * g + 1] = 0.0;
            v4tau4[5 * g + 2] = 0.0;
            v4tau4[5 * g + 3] = 0.0;
            v4tau4[5 * g + 4] = 0.0;
        }
    }
}

void
copyWeights(double* screenedWeights, const int32_t gridBlockPosition, const double* weights, const int32_t nGridPoints)
{
    for (int32_t g = 0; g < nGridPoints; g++)
    {
        screenedWeights[g] = weights[gridBlockPosition + g];
    }
}

}  // namespace gridscreen
