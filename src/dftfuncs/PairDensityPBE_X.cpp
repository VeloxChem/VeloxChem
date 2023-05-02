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

#include "PairDensityPBE_X.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftpbe_x {  // pdftpbe_x namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double one3 = 1.0 / 3.0;

    double two3 = 2.0 / 3.0;

    double four3 = 4.0 / 3.0;

    double five3 = 5.0 / 3.0;

    double eight3 = 8.0 / 3.0;

    double eleven3 = 11.0 / 3.0;

    double frg = -0.75 * std::pow(3.0 / mathconst::getPiValue(), one3);

    double R = 0.804;

    double Rpo = 1.0 + R;

    double mu = 0.066725 / 3.0 * std::pow(mathconst::getPiValue(), 2);

    double mus2r = mu * std::pow(1.0 / 12 * std::pow(6.0 / mathconst::getPiValue(), 2.0 / 3.0), 2) / R;

    for (int32_t g = 0; g < np; g++)
    {
        double density = rho[2 * g + 0];

        if (density < 1.0e-8)
        {
            exc[g] = 0.0;

            vrho[2 * g + 0] = 0.0;
            vrho[2 * g + 1] = 0.0;

            vsigma[3 * g + 0] = 0.0;
            vsigma[3 * g + 1] = 0.0;
            vsigma[3 * g + 2] = 0.0;

            continue;
        }

        double rho13 = std::pow(density, one3);

        double rho43 = std::pow(density, four3);

        double drho43 = four3 * std::pow(density, one3);

        double rho53 = std::pow(density, five3);

        double rho83 = std::pow(density, eight3);

        double rho2 = std::pow(density, 2.0);

        double pair_density = rho[2 * g + 1];

        double sig = sigma[3 * g + 0];

        double delta2 = -2.0 * pair_density;

        double delta = std::sqrt(std::fabs(-2.0 * pair_density));

        double zeta = delta / density;

        double zeta2 = zeta * zeta;  // This is really the absolute value

        // Li Manni's 2014 translation of gradients
        double gradR = 0.25 * sig * (1.0 + delta2 / (density * density));

        double gradI = 0.25 * sig * zeta;

        double f_zeta;

        double fg_zeta;

        double dExc_rho;

        double dExc_pi;

        double dExc_sigma;

        // Real case
        if (pair_density <= 0)
        {
            double fa = std::pow(1.0 + zeta, four3);

            double fb = std::pow(std::max(1.0 - zeta, 0.0), four3);

            f_zeta = (fa + fb);

            double f1_zeta = std::pow((1 + zeta), one3) - std::pow((1 - zeta), one3);

            double rhoa = 0.5 * (density + delta);

            double rhob = 0.5 * (density - delta);

            double grada2 = gradR + gradI;

            double gradb2 = gradR - gradI;

            double alpha = 1.0 + mus2r * grada2 / std::pow(rhoa, eight3);

            double beta = 1.0 + mus2r * gradb2 / std::pow(rhob, eight3);

            double Fxc_a = R / alpha;

            double Fxc_b = R;

            if (rhob > 1.0e-12)
            {
                Fxc_b = R / beta;
            }

            fg_zeta = fa * Fxc_a + fb * Fxc_b;

            //--- Derivatives
            // rho
            double dfa_rho     = -four3 * zeta / density * std::pow((1 + zeta), one3);
            double dfb_rho     = four3 * zeta / density * std::pow((1 - zeta), one3);
            double dF_rho      = -four3 * zeta / density * f1_zeta;
            double dgradR_rho  = -0.5 * sig * zeta2 / density;
            double dgradI_rho  = -0.25 * sig * delta / rho2;
            double dgrada2_rho = dgradR_rho + dgradI_rho;
            double dgradb2_rho = dgradR_rho - dgradI_rho;
            double dFxc_a_rho  = -R / std::pow(alpha, 2) * mus2r / std::pow(rhoa, eight3) * (dgrada2_rho - four3 / rhoa * grada2);
            double dFxc_b_rho  = -R / std::pow(beta, 2) * mus2r / std::pow(rhob, eight3) * (dgradb2_rho - four3 / rhob * gradb2);
            double dFg_rho     = (dfa_rho * Fxc_a + fa * dFxc_a_rho) + (dfb_rho * Fxc_b + fb * dFxc_b_rho);

            dExc_rho = 0.5 * frg * (Rpo * (drho43 * f_zeta + rho43 * dF_rho) - (drho43 * fg_zeta + rho43 * dFg_rho));

            // pi
            double dfa_pi     = -four3 / (density * delta) * std::pow((1 + zeta), one3);
            double dfb_pi     = four3 / (density * delta) * std::pow((1 - zeta), one3);
            double dF_pi      = dfa_pi + dfb_pi;
            double dgradR_pi  = -0.5 * sig / (density * density);
            double dgradI_pi  = -0.25 * sig / (density * delta);
            double dgrada2_pi = dgradR_pi + dgradI_pi;
            double dgradb2_pi = dgradR_pi - dgradI_pi;
            double dFxc_a_pi  = -R / std::pow(alpha, 2.0) * mus2r / std::pow(rhoa, eight3) * (dgrada2_pi + four3 * grada2 / rhoa / delta);
            double dFxc_b_pi  = -R / std::pow(beta, 2.0) * mus2r / std::pow(rhob, eight3) * (dgradb2_pi - four3 * gradb2 / rhob / delta);
            double dFg_pi     = (dfa_pi * Fxc_a + fa * dFxc_a_pi) + (dfb_pi * Fxc_b + fb * dFxc_b_pi);

            dExc_pi = 0.5 * frg * rho43 * (Rpo * dF_pi - dFg_pi);

            // sigma
            double dgradR_sigma  = 0.25 * (1 + std::pow(zeta, 2));
            double dgradI_sigma  = 0.25 * zeta;
            double dgrada2_sigma = dgradR_sigma + dgradI_sigma;
            double dgradb2_sigma = dgradR_sigma - dgradI_sigma;
            double dFxc_a_sigma  = -R / std::pow(alpha, 2.0) * mus2r / std::pow(rhoa, eight3) * dgrada2_sigma;
            double dFxc_b_sigma  = -R / std::pow(beta, 2.0) * mus2r / std::pow(rhob, eight3) * dgradb2_sigma;
            double dFg_sigma     = fa * dFxc_a_sigma + fb * dFxc_b_sigma;

            dExc_sigma = -0.5 * frg * rho43 * dFg_sigma;
        }
        // Imaginary case
        else
        {
            double r = std::sqrt(1.0 + zeta2);

            double theta = 4.0 / 3.0 * std::atan(zeta);

            f_zeta = 2.0 * std::pow(r, 4.0 / 3.0) * std::cos(theta);

            double f1_zeta = zeta2 * std::cos(theta) - zeta * std::sin(theta);

            double rho83 = std::pow(density, 8.0 / 3.0);

            double denom_R = rho83 * std::pow(0.5 * r, 8.0 / 3.0) * std::cos(2.0 * theta) + mus2r * gradR;

            double denom_I = rho83 * std::pow(0.5 * r, 8.0 / 3.0) * std::sin(2.0 * theta) + mus2r * gradI;

            double r_denom = std::sqrt(denom_R * denom_R + denom_I * denom_I);

            double theta_denom = std::atan2(denom_I, denom_R);

            double r_final = std::pow(0.5, 8.0 / 3.0) * std::pow(r, 4.0) / r_denom;

            double theta_final = 3.0 * theta - theta_denom;

            fg_zeta = 2.0 * R * rho83 * r_final * std::cos(theta_final);

            //--- Derivatives
            // rho
            double dr_rho     = -std::pow(zeta, 2.0) / (r * density);
            double dtheta_rho = -four3 * zeta / (density * std::pow(r, 2.0));
            double dF_rho     = -eight3 / (density * std::pow(r, two3)) * f1_zeta;
            double dgradR_rho = -0.5 * sig * zeta2 / density;
            double dgradI_rho = -0.25 * sig * delta / rho2;
            double dmR_rho    = eight3 * std::pow(0.5, eight3) * std::pow(density * r, five3) * (r + density * dr_rho) * std::cos(2 * theta) +
                             std::pow(density * r / 2.0, eight3) * (-2 * std::sin(2 * theta)) * dtheta_rho + mus2r * dgradR_rho;
            double dmI_rho = eight3 * std::pow(0.5, eight3) * std::pow(density * r, five3) * (r + density * dr_rho) * std::sin(2 * theta) +
                             std::pow(density * r / 2, eight3) * (2 * std::cos(2 * theta)) * dtheta_rho + mus2r * dgradI_rho;
            double drm_rho      = 1.0 / r_denom * (denom_R * dmR_rho + denom_I * dmI_rho);
            double dtheta_m_rho = 1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_rho / denom_I - denom_R / std::pow(denom_I, 2) * dmI_rho);
            double drf_rho = std::pow(0.5, eight3) * (4.0 * std::pow(r, 3.0) / r_denom * dr_rho - std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_rho);
            double dtheta_f_rho = 3.0 * dtheta_rho - dtheta_m_rho;
            double dFg_rho =
                2.0 * R *
                ((eight3 * rho53 * r_final + rho83 * drf_rho) * std::cos(theta_final) - rho83 * r_final * std::sin(theta_final) * dtheta_f_rho);

            dExc_rho = 0.5 * frg * (Rpo * (drho43 * f_zeta + rho43 * dF_rho) - (drho43 * fg_zeta + rho43 * dFg_rho));

            // pi
            double dr_pi     = -1.0 / (r * std::pow(density, 2.0));
            double dtheta_pi = -four3 / (std::pow(r, 2.0) * density * delta);
            double dF_pi     = eight3 * std::pow(r, one3) * std::cos(theta) * dr_pi - 2 * std::pow(r, four3) * std::sin(theta) * dtheta_pi;
            double dgradR_pi = -0.5 * sig / (density * density);
            double dgradI_pi = -0.25 * sig / (density * delta);
            double dmR_pi    = std::pow(density / 2.0, eight3) * (eight3 * std::pow(r, five3) * std::cos(2 * theta) * dr_pi -
                                                               2.0 * std::pow(r, eight3) * std::sin(2 * theta) * dtheta_pi) +
                            mus2r * dgradR_pi;
            double dmI_pi = std::pow(density / 2.0, eight3) * (eight3 * std::pow(r, five3) * std::sin(2 * theta) * dr_pi +
                                                               2.0 * std::pow(r, eight3) * std::cos(2 * theta) * dtheta_pi) +
                            mus2r * dgradI_pi;
            double drm_pi      = 1.0 / r_denom * (denom_R * dmR_pi + denom_I * dmI_pi);
            double dtheta_m_pi = 1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_pi / denom_I - denom_R / std::pow(denom_I, 2) * dmI_pi);
            double drf_pi      = std::pow(0.5, eight3) * (4 * std::pow(r, 3.0) / r_denom * dr_pi - std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_pi);
            double dtheta_f_pi = 3.0 * dtheta_pi - dtheta_m_pi;
            double dFg_pi      = 2.0 * R * rho83 * (drf_pi * std::cos(theta_final) - r_final * std::sin(theta_final) * dtheta_f_pi);

            dExc_pi = -0.5 * frg * rho43 * (Rpo * dF_pi - dFg_pi);

            // sigma
            double dgradR_sigma   = 0.25 * (1 + std::pow(zeta, 2));
            double dgradI_sigma   = 0.25 * zeta;
            double dmR_sigma      = mus2r * dgradR_sigma;
            double dmI_sigma      = mus2r * dgradI_sigma;
            double drm_sigma      = 1.0 / r_denom * (denom_R * dmR_sigma + denom_I * dmI_sigma);
            double dtheta_m_sigma = 1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_sigma / denom_I - denom_R / std::pow(denom_I, 2) * dmI_sigma);
            double drf_sigma      = std::pow(0.5, eight3) * (-std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_sigma);
            double dtheta_f_sigma = -dtheta_m_sigma;
            double dFg_sigma      = 2.0 * R * rho83 * (drf_sigma * std::cos(theta_final) - r_final * std::sin(theta_final) * dtheta_f_sigma);

            dExc_sigma = -0.5 * frg * rho43 * dFg_sigma;
        }

        exc[g] = 0.5 * frg * rho13 * (f_zeta * (1.0 + R) - fg_zeta);

        vrho[2 * g + 0] = dExc_rho;
        vrho[2 * g + 1] = dExc_pi;

        vsigma[3 * g + 0] = dExc_sigma;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_x
