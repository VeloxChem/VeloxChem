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
screenExcVxcForLDA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, double* exc, double* vrho)
{
    double densityThreshold = getDensityScreeningThreshold();

    auto ldafunc = xcFunctionalPointer->getFunctionalPointerToLdaComponent();

    const auto dim = &(ldafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho
        if (((std::fabs(rho[2 * g + 0]) <= densityThreshold) && (std::fabs(rho[2 * g + 1]) <= densityThreshold)))
        {
            exc[g] = 0.0;
        }

        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
        }
    }
}

void
screenVxcForLDA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, double* vrho)
{
    double densityThreshold = getDensityScreeningThreshold();

    auto ldafunc = xcFunctionalPointer->getFunctionalPointerToLdaComponent();

    const auto dim = &(ldafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
        }
    }
}

void
screenExcVxcForGGA(const CXCFunctional* xcFunctionalPointer,
                   const int32_t        npoints,
                   const double*        rho,
                   const double*        sigma,
                   double*              exc,
                   double*              vrho,
                   double*              vsigma)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    auto ggafunc = xcFunctionalPointer->getFunctionalPointerToGgaComponent();

    const auto dim = &(ggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
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
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vsigma - 1; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vsigma; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
        }
    }
}

void
screenVxcForGGA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, const double* sigma, double* vrho, double* vsigma)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    auto ggafunc = xcFunctionalPointer->getFunctionalPointerToGgaComponent();

    const auto dim = &(ggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vsigma - 1; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vsigma; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
        }
    }
}

void
screenExcVxcForMGGA(const CXCFunctional* xcFunctionalPointer,
                    const int32_t npoints,
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

    auto mggafunc = xcFunctionalPointer->getFunctionalPointerToMetaGgaComponent();

    const auto dim = &(mggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
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
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vsigma - 1; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vlapl - 1; ind++)
            {
                vlapl[dim->vlapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vtau - 1; ind++)
            {
                vtau[dim->vtau * g + ind] = 0.0;
            }
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vsigma; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vlapl; ind++)
            {
                vlapl[dim->vlapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vtau; ind++)
            {
                vtau[dim->vtau * g + ind] = 0.0;
            }
        }
    }
}

void
screenVxcForMGGA(const CXCFunctional* xcFunctionalPointer,
                 const int32_t npoints,
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

    auto mggafunc = xcFunctionalPointer->getFunctionalPointerToMetaGgaComponent();

    const auto dim = &(mggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            for (int ind = 0; ind < dim->vrho - 1; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vsigma - 1; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vlapl - 1; ind++)
            {
                vlapl[dim->vlapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->vtau - 1; ind++)
            {
                vtau[dim->vtau * g + ind] = 0.0;
            }
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            for (int ind = 1; ind < dim->vrho; ind++)
            {
                vrho[dim->vrho * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vsigma; ind++)
            {
                vsigma[dim->vsigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vlapl; ind++)
            {
                vlapl[dim->vlapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->vtau; ind++)
            {
                vtau[dim->vtau * g + ind] = 0.0;
            }
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
screenFxcForLDA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, double* v2rho2)
{
    double densityThreshold = getDensityScreeningThreshold();

    auto ldafunc = xcFunctionalPointer->getFunctionalPointerToLdaComponent();

    const auto dim = &(ldafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            for (int ind = 0; ind < dim->v2rho2 - 1; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            for (int ind = 1; ind < dim->v2rho2; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
        }
    }
}

void
screenFxcForGGA(const CXCFunctional* xcFunctionalPointer,
                const int32_t        npoints,
                const double*        rho,
                const double*        sigma,
                double*              v2rho2,
                double*              v2rhosigma,
                double*              v2sigma2)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    auto ggafunc = xcFunctionalPointer->getFunctionalPointerToGgaComponent();

    const auto dim = &(ggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            for (int ind = 0; ind < dim->v2rho2 - 1; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2rhosigma - 1; ind++)
            {
                v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2sigma2 - 1; ind++)
            {
                v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
            }
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            for (int ind = 1; ind < dim->v2rho2; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2rhosigma; ind++)
            {
                v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2sigma2; ind++)
            {
                v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
            }
        }
    }
}

void
screenFxcForMGGA(const CXCFunctional* xcFunctionalPointer,
                 const int32_t        npoints,
                 const double*        rho,
                 const double*        sigma,
                 const double*        lapl,
                 const double*        tau,
                 double*              v2rho2,
                 double*              v2rhosigma,
                 double*              v2rholapl,
                 double*              v2rhotau,
                 double*              v2sigma2,
                 double*              v2sigmalapl,
                 double*              v2sigmatau,
                 double*              v2lapl2,
                 double*              v2lapltau,
                 double*              v2tau2)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    auto mggafunc = xcFunctionalPointer->getFunctionalPointerToMetaGgaComponent();

    const auto dim = &(mggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            for (int ind = 0; ind < dim->v2rho2 - 1; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2rhosigma - 1; ind++)
            {
                v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2rholapl - 1; ind++)
            {
                v2rholapl[dim->v2rholapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2rhotau - 1; ind++)
            {
                v2rhotau[dim->v2rhotau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2sigma2 - 1; ind++)
            {
                v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2sigmalapl - 1; ind++)
            {
                v2sigmalapl[dim->v2sigmalapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2sigmatau - 1; ind++)
            {
                v2sigmatau[dim->v2sigmatau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2lapl2 - 1; ind++)
            {
                v2lapl2[dim->v2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2lapltau - 1; ind++)
            {
                v2lapltau[dim->v2lapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v2tau2 - 1; ind++)
            {
                v2tau2[dim->v2tau2 * g + ind] = 0.0;
            }
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            for (int ind = 1; ind < dim->v2rho2; ind++)
            {
                v2rho2[dim->v2rho2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2rhosigma; ind++)
            {
                v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2rholapl; ind++)
            {
                v2rholapl[dim->v2rholapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2rhotau; ind++)
            {
                v2rhotau[dim->v2rhotau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2sigma2; ind++)
            {
                v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2sigmalapl; ind++)
            {
                v2sigmalapl[dim->v2sigmalapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2sigmatau; ind++)
            {
                v2sigmatau[dim->v2sigmatau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2lapl2; ind++)
            {
                v2lapl2[dim->v2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2lapltau; ind++)
            {
                v2lapltau[dim->v2lapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v2tau2; ind++)
            {
                v2tau2[dim->v2tau2 * g + ind] = 0.0;
            }
        }
    }
}

void
screenKxcForLDA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, double* v3rho3)
{
    double densityThreshold = getDensityScreeningThreshold();

    auto ldafunc = xcFunctionalPointer->getFunctionalPointerToLdaComponent();

    const auto dim = &(ldafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            for (int ind = 0; ind < dim->v3rho3 - 1; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            for (int ind = 1; ind < dim->v3rho3; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
        }
    }
}

void
screenKxcForGGA(const CXCFunctional* xcFunctionalPointer,
                const int32_t        npoints,
                const double*        rho,
                const double*        sigma,
                double*              v3rho3,
                double*              v3rho2sigma,
                double*              v3rhosigma2,
                double*              v3sigma3)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    auto ggafunc = xcFunctionalPointer->getFunctionalPointerToGgaComponent();

    const auto dim = &(ggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            for (int ind = 0; ind < dim->v3rho3 - 1; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rho2sigma - 1; ind++)
            {
                v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rhosigma2 - 1; ind++)
            {
                v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigma3 - 1; ind++)
            {
                v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
            }
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            for (int ind = 1; ind < dim->v3rho3; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rho2sigma; ind++)
            {
                v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rhosigma2; ind++)
            {
                v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigma3; ind++)
            {
                v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
            }
        }
    }
}

void
screenKxcForMGGA(const CXCFunctional* xcFunctionalPointer,
                 const int32_t        npoints,
                 const double*        rho,
                 const double*        sigma,
                 const double*        lapl,
                 const double*        tau,
                 double*              v3rho3,
                 double*              v3rho2sigma,
                 double*              v3rho2lapl,
                 double*              v3rho2tau,
                 double*              v3rhosigma2,
                 double*              v3rhosigmalapl,
                 double*              v3rhosigmatau,
                 double*              v3rholapl2,
                 double*              v3rholapltau,
                 double*              v3rhotau2,
                 double*              v3sigma3,
                 double*              v3sigma2lapl,
                 double*              v3sigma2tau,
                 double*              v3sigmalapl2,
                 double*              v3sigmalapltau,
                 double*              v3sigmatau2,
                 double*              v3lapl3,
                 double*              v3lapl2tau,
                 double*              v3lapltau2,
                 double*              v3tau3)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    auto mggafunc = xcFunctionalPointer->getFunctionalPointerToMetaGgaComponent();

    const auto dim = &(mggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            for (int ind = 0; ind < dim->v3rho3 - 1; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rho2sigma - 1; ind++)
            {
                v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rho2lapl - 1; ind++)
            {
                v3rho2lapl[dim->v3rho2lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rho2tau - 1; ind++)
            {
                v3rho2tau[dim->v3rho2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rhosigma2 - 1; ind++)
            {
                v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rhosigmalapl - 1; ind++)
            {
                v3rhosigmalapl[dim->v3rhosigmalapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rhosigmatau - 1; ind++)
            {
                v3rhosigmatau[dim->v3rhosigmatau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rholapl2 - 1; ind++)
            {
                v3rholapl2[dim->v3rholapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rholapltau - 1; ind++)
            {
                v3rholapltau[dim->v3rholapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3rhotau2 - 1; ind++)
            {
                v3rhotau2[dim->v3rhotau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigma3 - 1; ind++)
            {
                v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigma2lapl - 1; ind++)
            {
                v3sigma2lapl[dim->v3sigma2lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigma2tau - 1; ind++)
            {
                v3sigma2tau[dim->v3sigma2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigmalapl2 - 1; ind++)
            {
                v3sigmalapl2[dim->v3sigmalapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigmalapltau - 1; ind++)
            {
                v3sigmalapltau[dim->v3sigmalapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3sigmatau2 - 1; ind++)
            {
                v3sigmatau2[dim->v3sigmatau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3lapl3 - 1; ind++)
            {
                v3lapl3[dim->v3lapl3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3lapl2tau - 1; ind++)
            {
                v3lapl2tau[dim->v3lapl2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3lapltau2 - 1; ind++)
            {
                v3lapltau2[dim->v3lapltau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v3tau3 - 1; ind++)
            {
                v3tau3[dim->v3tau3 * g + ind] = 0.0;
            }
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            for (int ind = 1; ind < dim->v3rho3; ind++)
            {
                v3rho3[dim->v3rho3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rho2sigma; ind++)
            {
                v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rho2lapl; ind++)
            {
                v3rho2lapl[dim->v3rho2lapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rho2tau; ind++)
            {
                v3rho2tau[dim->v3rho2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rhosigma2; ind++)
            {
                v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rhosigmalapl; ind++)
            {
                v3rhosigmalapl[dim->v3rhosigmalapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rhosigmatau; ind++)
            {
                v3rhosigmatau[dim->v3rhosigmatau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rholapl2; ind++)
            {
                v3rholapl2[dim->v3rholapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rholapltau; ind++)
            {
                v3rholapltau[dim->v3rholapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3rhotau2; ind++)
            {
                v3rhotau2[dim->v3rhotau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigma3; ind++)
            {
                v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigma2lapl; ind++)
            {
                v3sigma2lapl[dim->v3sigma2lapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigma2tau; ind++)
            {
                v3sigma2tau[dim->v3sigma2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigmalapl2; ind++)
            {
                v3sigmalapl2[dim->v3sigmalapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigmalapltau; ind++)
            {
                v3sigmalapltau[dim->v3sigmalapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3sigmatau2; ind++)
            {
                v3sigmatau2[dim->v3sigmatau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3lapl3; ind++)
            {
                v3lapl3[dim->v3lapl3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3lapl2tau; ind++)
            {
                v3lapl2tau[dim->v3lapl2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3lapltau2; ind++)
            {
                v3lapltau2[dim->v3lapltau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v3tau3; ind++)
            {
                v3tau3[dim->v3tau3 * g + ind] = 0.0;
            }
        }
    }
}

void
screenLxcForLDA(const CXCFunctional* xcFunctionalPointer, const int32_t npoints, const double* rho, double* v4rho4)
{
    double densityThreshold = getDensityScreeningThreshold();

    auto ldafunc = xcFunctionalPointer->getFunctionalPointerToLdaComponent();

    const auto dim = &(ldafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a
        if (std::fabs(rho[2 * g + 0]) <= densityThreshold)
        {
            for (int ind = 0; ind < dim->v4rho4 - 1; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
        }

        // rho_b
        if (std::fabs(rho[2 * g + 1]) <= densityThreshold)
        {
            for (int ind = 1; ind < dim->v4rho4; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
        }
    }
}

void
screenLxcForGGA(const CXCFunctional* xcFunctionalPointer,
                const int32_t        npoints,
                const double*        rho,
                const double*        sigma,
                double*              v4rho4,
                double*              v4rho3sigma,
                double*              v4rho2sigma2,
                double*              v4rhosigma3,
                double*              v4sigma4)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    auto ggafunc = xcFunctionalPointer->getFunctionalPointerToGgaComponent();

    const auto dim = &(ggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a and sigma_aa
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold))
        {
            for (int ind = 0; ind < dim->v4rho4 - 1; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3sigma - 1; ind++)
            {
                v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigma2 - 1; ind++)
            {
                v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma3 - 1; ind++)
            {
                v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma4 - 1; ind++)
            {
                v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
            }
        }

        // rho_b and sigma_bb
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold))
        {
            for (int ind = 1; ind < dim->v4rho4; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho3sigma; ind++)
            {
                v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2sigma2; ind++)
            {
                v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigma3; ind++)
            {
                v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma4; ind++)
            {
                v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
            }
        }
    }
}

void
screenLxcForMGGA(const CXCFunctional* xcFunctionalPointer,
                 const int32_t        npoints,
                 const double*        rho,
                 const double*        sigma,
                 const double*        lapl,
                 const double*        tau,
                 double*              v4rho4,
                 double*              v4rho3sigma,
                 double*              v4rho3lapl,
                 double*              v4rho3tau,
                 double*              v4rho2sigma2,
                 double*              v4rho2sigmalapl,
                 double*              v4rho2sigmatau,
                 double*              v4rho2lapl2,
                 double*              v4rho2lapltau,
                 double*              v4rho2tau2,
                 double*              v4rhosigma3,
                 double*              v4rhosigma2lapl,
                 double*              v4rhosigma2tau,
                 double*              v4rhosigmalapl2,
                 double*              v4rhosigmalapltau,
                 double*              v4rhosigmatau2,
                 double*              v4rholapl3,
                 double*              v4rholapl2tau,
                 double*              v4rholapltau2,
                 double*              v4rhotau3,
                 double*              v4sigma4,
                 double*              v4sigma3lapl,
                 double*              v4sigma3tau,
                 double*              v4sigma2lapl2,
                 double*              v4sigma2lapltau,
                 double*              v4sigma2tau2,
                 double*              v4sigmalapl3,
                 double*              v4sigmalapl2tau,
                 double*              v4sigmalapltau2,
                 double*              v4sigmatau3,
                 double*              v4lapl4,
                 double*              v4lapl3tau,
                 double*              v4lapl2tau2,
                 double*              v4lapltau3,
                 double*              v4tau4)
{
    double densityThreshold = getDensityScreeningThreshold();

    double sigmaThreshold = getSigmaScreeningThreshold(densityThreshold);

    double tauThreshold = getTauScreeningThreshold();

    auto mggafunc = xcFunctionalPointer->getFunctionalPointerToMetaGgaComponent();

    const auto dim = &(mggafunc->dim);

    for (int g = 0; g < static_cast<int>(npoints); g++)
    {
        // rho_a, sigma_aa and tau_a
        if ((std::fabs(rho[2 * g + 0]) <= densityThreshold) || (std::fabs(sigma[3 * g + 0]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 0]) <= tauThreshold))
        {
            for (int ind = 0; ind < dim->v4rho4 - 1; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3sigma - 1; ind++)
            {
                v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3lapl - 1; ind++)
            {
                v4rho3lapl[dim->v4rho3lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3tau - 1; ind++)
            {
                v4rho3tau[dim->v4rho3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigma2 - 1; ind++)
            {
                v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigmalapl - 1; ind++)
            {
                v4rho2sigmalapl[dim->v4rho2sigmalapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigmatau - 1; ind++)
            {
                v4rho2sigmatau[dim->v4rho2sigmatau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2lapl2 - 1; ind++)
            {
                v4rho2lapl2[dim->v4rho2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2lapltau - 1; ind++)
            {
                v4rho2lapltau[dim->v4rho2lapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2tau2 - 1; ind++)
            {
                v4rho2tau2[dim->v4rho2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma3 - 1; ind++)
            {
                v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma2lapl - 1; ind++)
            {
                v4rhosigma2lapl[dim->v4rhosigma2lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma2tau - 1; ind++)
            {
                v4rhosigma2tau[dim->v4rhosigma2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmalapl2 - 1; ind++)
            {
                v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmalapltau - 1; ind++)
            {
                v4rhosigmalapltau[dim->v4rhosigmalapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmatau2 - 1; ind++)
            {
                v4rhosigmatau2[dim->v4rhosigmatau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapl3 - 1; ind++)
            {
                v4rholapl3[dim->v4rholapl3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapl2tau - 1; ind++)
            {
                v4rholapl2tau[dim->v4rholapl2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapltau2 - 1; ind++)
            {
                v4rholapltau2[dim->v4rholapltau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhotau3 - 1; ind++)
            {
                v4rhotau3[dim->v4rhotau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma4 - 1; ind++)
            {
                v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma3lapl - 1; ind++)
            {
                v4sigma3lapl[dim->v4sigma3lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma3tau - 1; ind++)
            {
                v4sigma3tau[dim->v4sigma3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2lapl2 - 1; ind++)
            {
                v4sigma2lapl2[dim->v4sigma2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2lapltau - 1; ind++)
            {
                v4sigma2lapltau[dim->v4sigma2lapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2tau2 - 1; ind++)
            {
                v4sigma2tau2[dim->v4sigma2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapl3 - 1; ind++)
            {
                v4sigmalapl3[dim->v4sigmalapl3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapl2tau - 1; ind++)
            {
                v4sigmalapl2tau[dim->v4sigmalapl2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapltau2 - 1; ind++)
            {
                v4sigmalapltau2[dim->v4sigmalapltau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmatau3 - 1; ind++)
            {
                v4sigmatau3[dim->v4sigmatau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl4 - 1; ind++)
            {
                v4lapl4[dim->v4lapl4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl3tau - 1; ind++)
            {
                v4lapl3tau[dim->v4lapl3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl2tau2 - 1; ind++)
            {
                v4lapl2tau2[dim->v4lapl2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapltau3 - 1; ind++)
            {
                v4lapltau3[dim->v4lapltau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4tau4 - 1; ind++)
            {
                v4tau4[dim->v4tau4 * g + ind] = 0.0;
            }
        }

        // rho_b, sigma_bb and tau_b
        if ((std::fabs(rho[2 * g + 1]) <= densityThreshold) || (std::fabs(sigma[3 * g + 2]) <= sigmaThreshold) ||
            (std::fabs(tau[2 * g + 1]) <= tauThreshold))
        {
            for (int ind = 1; ind < dim->v4rho4; ind++)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho3sigma; ind++)
            {
                v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho3lapl; ind++)
            {
                v4rho3lapl[dim->v4rho3lapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho3tau; ind++)
            {
                v4rho3tau[dim->v4rho3tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2sigma2; ind++)
            {
                v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2sigmalapl; ind++)
            {
                v4rho2sigmalapl[dim->v4rho2sigmalapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2sigmatau; ind++)
            {
                v4rho2sigmatau[dim->v4rho2sigmatau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2lapl2; ind++)
            {
                v4rho2lapl2[dim->v4rho2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2lapltau; ind++)
            {
                v4rho2lapltau[dim->v4rho2lapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rho2tau2; ind++)
            {
                v4rho2tau2[dim->v4rho2tau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigma3; ind++)
            {
                v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigma2lapl; ind++)
            {
                v4rhosigma2lapl[dim->v4rhosigma2lapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigma2tau; ind++)
            {
                v4rhosigma2tau[dim->v4rhosigma2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigmalapl2; ind++)
            {
                v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigmalapltau; ind++)
            {
                v4rhosigmalapltau[dim->v4rhosigmalapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhosigmatau2; ind++)
            {
                v4rhosigmatau2[dim->v4rhosigmatau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rholapl3; ind++)
            {
                v4rholapl3[dim->v4rholapl3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rholapl2tau; ind++)
            {
                v4rholapl2tau[dim->v4rholapl2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rholapltau2; ind++)
            {
                v4rholapltau2[dim->v4rholapltau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4rhotau3; ind++)
            {
                v4rhotau3[dim->v4rhotau3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma4; ind++)
            {
                v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma3lapl; ind++)
            {
                v4sigma3lapl[dim->v4sigma3lapl * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma3tau; ind++)
            {
                v4sigma3tau[dim->v4sigma3tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma2lapl2; ind++)
            {
                v4sigma2lapl2[dim->v4sigma2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma2lapltau; ind++)
            {
                v4sigma2lapltau[dim->v4sigma2lapltau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigma2tau2; ind++)
            {
                v4sigma2tau2[dim->v4sigma2tau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigmalapl3; ind++)
            {
                v4sigmalapl3[dim->v4sigmalapl3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigmalapl2tau; ind++)
            {
                v4sigmalapl2tau[dim->v4sigmalapl2tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigmalapltau2; ind++)
            {
                v4sigmalapltau2[dim->v4sigmalapltau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4sigmatau3; ind++)
            {
                v4sigmatau3[dim->v4sigmatau3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4lapl4; ind++)
            {
                v4lapl4[dim->v4lapl4 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4lapl3tau; ind++)
            {
                v4lapl3tau[dim->v4lapl3tau * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4lapl2tau2; ind++)
            {
                v4lapl2tau2[dim->v4lapl2tau2 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4lapltau3; ind++)
            {
                v4lapltau3[dim->v4lapltau3 * g + ind] = 0.0;
            }
            for (int ind = 1; ind < dim->v4tau4; ind++)
            {
                v4tau4[dim->v4tau4 * g + ind] = 0.0;
            }
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
