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

#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftpbe_x {  // pdftpbe_x namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double onethird = 1.0 / 3.0;

    double fourthird = 4.0 / 3.0;

    double frg = - 0.75 * std::pow(3.0 / mathconst::getPiValue(), onethird);

    double R = 0.804;

    double mu = 0.066725 / 3.0 *  std::pow(mathconst::getPiValue(), 2);

    double mus2r = mu * pow(1.0 / 12 * std::pow( 6.0 /mathconst::getPiValue(), 2.0/ 3.0 ), 2) / R;

    //Translation parameters

    double transR0 = 0.01;

    double transgrad_A = 156250000;

    double transgrad_B = -56250;

    double transgrad_C = 14.0625;

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

        double pair_density = rho[2 * g + 1];

        double delta2 = -2.0 * pair_density;

        double rho13 = std::pow(density, onethird);

        double delta;

        if (pair_density <= 0)
        {
            delta = sqrt(-2.0 * pair_density);
        }
        else
        {
            delta = sqrt(2.0 * pair_density);
        }

        double zeta = delta / density;

        double zeta2 = zeta*zeta; //This is really the absolute value

        double invdelta;

        if (zeta2 < transR0)
        {
           invdelta = (transgrad_A * std::pow(zeta2, 4.0) + transgrad_B * std::pow(zeta2, 2.0) + transgrad_C) / density;
        }
        else
        {
            invdelta = 1.0 / delta;
        }

        double gradR = 0.25 * (sigma[3 * g + 0] + sigma[3 * g + 2] * invdelta * invdelta);

        double gradI = 0.5 * sigma[3 * g + 1] * invdelta;

        double f_zeta;

        double fg_zeta;

        // Real case
        if (pair_density <= 0)
        {
            double fa = pow(1.0 + zeta, fourthird);

            double fb = pow(std::max(1.0 - zeta,0.0), fourthird);

            f_zeta = (fa + fb) * (1.0 + R);

            double rhoa = 0.5*(density + delta);

            double rhob = 0.5*(density - delta);

            double grada2 = gradR - gradI;

            double gradb2 = gradR + gradI;

            double Fxc_a = R/(1.0 + mus2r * grada2/pow(rhoa,2.666666666666667));

            double Fxc_b = 0.0;

            if (rhob > 1.0e-12)
            {
                Fxc_b = R/(1.0 + mus2r * gradb2/pow(rhob,2.666666666666667));
            }

            fg_zeta = fa * Fxc_a + fb * Fxc_b;
        }
        // Imaginary case
        else
        {
            double r = sqrt(1.0 + zeta2);

            double theta = 4.0 / 3.0 * std::atan(zeta);

            f_zeta = 2.0 * std::pow(r, 4.0 / 3.0) * std::cos(theta) * (1.0 + R);

            double rho83 = std::pow (density, 8.0 / 3.0);

            double denom_R = rho83 * std::pow(0.5 * r, 8.0 / 3.0) * std::cos(2.0 * theta) + mus2r * gradR;

            double denom_I = rho83 * std::pow(0.5 * r, 8.0 / 3.0) * std::sin(2.0 * theta) + mus2r * gradI;

            double r_denom = sqrt( denom_R * denom_R + denom_I * denom_I);

            double theta_denom = std::atan2(denom_I, denom_R);

            double r_final = std::pow(0.5, 8.0 / 3.0) * std::pow(r, 4.0) / r_denom;

            double theta_final = 2.0 * theta - theta_denom;

            fg_zeta = 2.0 * R * rho83 * r_final * std::cos(theta_final);
        }

        exc[g] = 0.5 * frg * rho13 * (f_zeta - fg_zeta);

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_x
