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
    double f13 = 1.0 / 3.0;

    double f23 = 2.0 / 3.0;

    double f43 = 4.0 / 3.0;

    double f53 = 5.0 / 3.0;

    double f83 = 8.0 / 3.0;

    double frg = -0.75 * std::pow(3.0 / mathconst::getPiValue(), f13);

    double R = 0.804;

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

        double rho13 = std::pow(density, f13);

        double rho43 = std::pow(density, f43);

        double rho53 = std::pow(density, f53);

        double rho2 = std::pow(density, 2.0);

        double pair_density = rho[2 * g + 1];

        double sig = sigma[3 * g + 0];

        double delta2 = -2.0 * pair_density;

        double delta = std::sqrt(std::fabs(-2.0 * pair_density));

        double zeta = delta / density;

        double zeta2 = zeta * zeta;  // This is really the absolute value

        // Li Manni's 2014 translation of gradients
        double gradR = 0.25 * sig * (1.0 + delta2 / rho2);

        double gradI = 0.5 * sig * zeta;

        double dgradR_rho  = -0.5 * sig * zeta2 / density;

        double dgradI_rho = -0.5 * sig * delta / rho2;

        double f_zeta = 0.0;

        double fg_zeta = 0.0;

        double dF_rho = 0.0;

        double dFg_rho = 0.0;

        double dF_pi = 0.0;

        double dFg_pi = 0.0;

        double dFg_sigma = 0.0;


        // Real case
        if (pair_density <= 0)
        {
            double fa = std::pow(1.0 + zeta, f43);

            double fb = std::pow(1.0 - zeta, f43);

            f_zeta = (fa + fb);

            double f1_zeta = f43 * (std::pow((1 + zeta), f13) - std::pow((1 - zeta), f13));

            double rhoa = 0.5 * (density + delta);

            double rhob = 0.5 * (density - delta);

            double grada2 = gradR + gradI;

            double gradb2 = gradR - gradI;

            double alpha = 1.0 + mus2r * grada2 / std::pow(rhoa, f83);

            double beta = 1.0;

            if (rhob > 1.0e-16)
            {
                beta = 1.0 + mus2r * gradb2 / std::pow(rhob, f83);
            }

            double Fxc_a = R / alpha;

            double Fxc_b = R / beta;

            fg_zeta = fa * Fxc_a + fb * Fxc_b;

            // Derivatives: rho

            double dfa_rho     = -f43 * zeta / density * std::pow((1 + zeta), f13);

            double dfb_rho     = f43 * zeta / density * std::pow((1 - zeta), f13);

            dF_rho      = - zeta / density * f1_zeta;

            double dFxc_a_pref = -R / std::pow(alpha, 2.0) * mus2r / std::pow(rhoa, f83);

            double dFxc_a_rho  = dFxc_a_pref * (dgradR_rho + dgradI_rho - f43 / rhoa * grada2);

            double dFxc_b_pref = 0.0;

            double dFxc_b_rho  = 0.0;

            if (rhob > 1.0e-16)
            {
                dFxc_b_pref = -R / std::pow(beta, 2.0) * mus2r / std::pow(rhob, f83);

                dFxc_b_rho  = dFxc_b_pref * (dgradR_rho - dgradI_rho - f43 / rhob * gradb2);
            }

            dFg_rho     = (dfa_rho * Fxc_a + fa * dFxc_a_rho) + (dfb_rho * Fxc_b + fb * dFxc_b_rho);

            // pi

            if (pair_density < -1.0e-12)
            {
                double dfa_pi     = -f43 / (density * delta) * std::pow((1 + zeta), f13);

                double dfb_pi     = f43 / (density * delta) * std::pow((1 - zeta), f13);

                dF_pi      = dfa_pi + dfb_pi;

                double dgradR_pi  = -0.5 * sig / rho2;

                double dgradI_pi  = -0.5 * sig / (density * delta);

                double dFxc_a_pi  = dFxc_a_pref * (dgradR_pi + dgradI_pi + f43 * grada2 / rhoa / delta);

                double dFxc_b_pi = 0.0;

                if (rhob > 1.0e-16)
                {
                    dFxc_b_pi  = dFxc_b_pref * (dgradR_pi - dgradI_pi - f43 * gradb2 / rhob / delta);
                }

                // 2* f43 / rhoa  * 0.25 * sig /rho

                dFg_pi     = (dfa_pi * Fxc_a + fa * dFxc_a_pi) + (dfb_pi * Fxc_b + fb * dFxc_b_pi);
            }
            else
            {
                dF_pi      = -f43 * f23 / rho2;

                double rho83 = std::pow(0.5 * density, f83);

                double A = 1.0 + mus2r * 0.25 * sig / rho83;

                dFg_pi = 8.0 / 9.0 * R / A - R * mus2r * sig / std::pow(A,2) / rho83 + 25.0/36.0 * R * std::pow(mus2r * sig/rho83, 2) / std::pow(A,3);

                dFg_pi = - dFg_pi / rho2;
            }

            // sigma

            double dgradR_sigma  = 0.25 * (1.0 + zeta2);

            double dgradI_sigma  = 0.5 * zeta;

            double dFxc_a_sigma  =  dFxc_a_pref * (dgradR_sigma + dgradI_sigma);

            double dFxc_b_sigma  =  dFxc_b_pref * (dgradR_sigma - dgradI_sigma);

            dFg_sigma     = fa * dFxc_a_sigma + fb * dFxc_b_sigma;

        }
        // Imaginary case
        else
        {
            double r = std::sqrt(1.0 + zeta2);

            double theta = 4.0 / 3.0 * std::atan(zeta);

            f_zeta = 2.0 * std::pow(r, 4.0 / 3.0) * std::cos(theta);

            double f1_zeta = f83 * (zeta * std::cos(theta) - std::sin(theta)) / std::pow(r, f23) ;

            double rho83 = std::pow(density, 8.0 / 3.0);

            double denom_pref = rho83 * std::pow(0.5 * r, 8.0 / 3.0);

            double denom_R = denom_pref * std::cos(2.0 * theta) + mus2r * gradR;

            double denom_I = denom_pref * std::sin(2.0 * theta) + mus2r * gradI;

            double r_denom = std::sqrt(denom_R * denom_R + denom_I * denom_I);

            double theta_denom = std::atan2(denom_I, denom_R);

            double r_final = std::pow(0.5, 8.0 / 3.0) * std::pow(r, 4.0) / r_denom;

            double theta_final = 3.0 * theta - theta_denom;

            fg_zeta = 2.0 * R * rho83 * r_final * std::cos(theta_final);

            // Derivatives: rho
            dF_rho     = - f1_zeta * zeta / density;

            double dr_rho     = - std::pow(zeta, 2.0) / (r * density);

            double dtheta_rho = -f43 * zeta / (density * std::pow(r, 2.0));

            double ddm_pref = 0.5 * std::pow(0.5 * density * r , f53);

            double dmR_rho = f83 * ddm_pref * (r + density * dr_rho) * std::cos(2 * theta) + denom_pref * (-2 * std::sin(2 * theta)) * dtheta_rho + mus2r * dgradR_rho;

            double dmI_rho = f83 * ddm_pref * (r + density * dr_rho) * std::sin(2 * theta) + denom_pref * (2 * std::cos(2 * theta)) * dtheta_rho + mus2r * dgradI_rho;

            double drm_rho      = 1.0 / r_denom * (denom_R * dmR_rho + denom_I * dmI_rho);

            double dtheta_m_rho = 1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_rho / denom_I - denom_R / std::pow(denom_I, 2) * dmI_rho);

            double drf_rho = std::pow(0.5, f83) * (4.0 * std::pow(r, 3.0) / r_denom * dr_rho - std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_rho);

            double dtheta_f_rho = 3.0 * dtheta_rho - dtheta_m_rho;

            dFg_rho = 2.0 * R *
                ((f83 * rho53 * r_final + rho83 * drf_rho) * std::cos(theta_final) - rho83 * r_final * std::sin(theta_final) * dtheta_f_rho);

            // pi
            if (pair_density > 1.0e-12)
            {
                double dr_pi     = 1.0 / (r * std::pow(density, 2.0));

                double dtheta_pi = f43 / (std::pow(r, 2.0) * density * delta);

                dF_pi     = f83 * std::pow(r, f13) * std::cos(theta) * dr_pi - 2 * std::pow(r, f43) * std::sin(theta) * dtheta_pi;

                double dgradR_pi  = 0.5 * sig / rho2;

                double dgradI_pi = 0.5 * sig / (density * delta);

                double dmR_pi    = ddm_pref * density * (f83 * std::cos(2 * theta) * dr_pi) - 2.0 * denom_pref * std::sin(2 * theta) * dtheta_pi + mus2r * dgradR_pi;

                double dmI_pi    = ddm_pref * density * (f83 * std::sin(2 * theta) * dr_pi) + 2.0 * denom_pref * std::cos(2 * theta) * dtheta_pi + mus2r * dgradI_pi;

                double drm_pi      = 1.0 / r_denom * (denom_R * dmR_pi + denom_I * dmI_pi);

                double dtheta_m_pi = -1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_pi / denom_I - denom_R / std::pow(denom_I, 2) * dmI_pi);

                double drf_pi      = std::pow(0.5, f83) * (4 * std::pow(r, 3.0) / r_denom * dr_pi - std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_pi);

                double dtheta_f_pi = 3.0 * dtheta_pi - dtheta_m_pi;

                dFg_pi      = 2.0 * R * rho83 * (drf_pi * std::cos(theta_final) - r_final * std::sin(theta_final) * dtheta_f_pi);
            }
            else
            {
                dF_pi      = -f43 * f23 / rho2;

                double rho83_2 = std::pow(0.5 * density, f83);

                double A = 1.0 + mus2r * 0.25 * sig / rho83_2;

                dFg_pi = 8.0 / 9.0 * R / A - R * mus2r * sig / std::pow(A,2) / rho83_2 + 25.0/36.0 * R * std::pow(mus2r * sig/rho83_2, 2) / std::pow(A,3);

                dFg_pi = - dFg_pi / rho2;
            }

            // sigma

            double dgradR_sigma   = 0.25 * (1 + std::pow(zeta, 2));

            double dgradI_sigma   = 0.5 * zeta;

            double dmR_sigma      = mus2r * dgradR_sigma;

            double dmI_sigma      = mus2r * dgradI_sigma;

            double drm_sigma      = 1.0 / r_denom * (denom_R * dmR_sigma + denom_I * dmI_sigma);

            double dtheta_m_sigma = 1.0 / (std::pow(denom_R / denom_I, 2) + 1) * (dmR_sigma / denom_I - denom_R / std::pow(denom_I, 2) * dmI_sigma);

            double drf_sigma      = std::pow(0.5, f83) * (-std::pow(r, 4.0) / std::pow(r_denom, 2) * drm_sigma);

            double dtheta_f_sigma = -dtheta_m_sigma;

            dFg_sigma      = 2.0 * R * rho83 * (drf_sigma * std::cos(theta_final) - r_final * std::sin(theta_final) * dtheta_f_sigma);
        }

        exc[g] = 0.5 * frg * rho13 * (f_zeta * (1.0 + R) - fg_zeta);

        vrho[2 * g + 0] = 0.5 * frg * rho43 * (dF_rho * (1.0 + R) - dFg_rho) + f43 * exc[g];
        vrho[2 * g + 1] = 0.5 * frg * rho43 * (dF_pi * (1.0 + R) - dFg_pi);

        vsigma[3 * g + 0] = -0.5 * frg * rho43 * dFg_sigma;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_x
