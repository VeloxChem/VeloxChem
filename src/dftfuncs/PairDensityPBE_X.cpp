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

        double delta = 0.0;

        if (pair_density < 0)
        {
            delta = sqrt(-2.0 * pair_density);
        }

        double zeta = delta / density;

        double rhoa = 0.5*(density + delta);

        double rhob = 0.5*(density - delta);

        double Fxc_a = 0.0 ; 

        double Fxc_b = 0.0 ; 

        if (delta != 0.0)
        {
            double grada2 = 0.25 * (sigma[3 * g + 0] - 2.0 * sigma[3 * g + 1]/delta + sigma[3 * g + 2]/delta2);

            double gradb2 = 0.25 * (sigma[3 * g + 0] + 2.0 * sigma[3 * g + 1]/delta + sigma[3 * g + 2]/delta2);


            Fxc_a = R/(1.0 + mus2r * grada2/pow(rhoa,2.666666666666667)); 

            if (rhob > 1.0e-12)
            {
                Fxc_b = R/(1.0 + mus2r * gradb2/pow(rhob,2.666666666666667));
            }
        }

        double fa = pow(1.0 + zeta, fourthird)*(1.0 + R - Fxc_a);

        double fb = pow(std::max(1.0 - zeta,0.0), fourthird)*(1.0 + R - Fxc_b);

        exc[g] = 0.5 * frg * rho13 * (fa + fb);

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_x
