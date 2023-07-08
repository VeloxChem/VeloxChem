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

#include "PairDensityLYP.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftlyp {  // pdftlyp namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
        
    double cf = 0.3 * std::pow(3 * mathconst::getPiValue() * mathconst::getPiValue(), 2.0 / 3.0);

    double f43 = 4.0/3.0;

    double f83 = 8.0/3.0;

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

        // Spin-polarization dependence
        double pair_density = rho[2 * g + 1];

        double density_gradient = sigma[3 * g + 0];

        double rho2 = density * density;

        double rhom13 = std::pow(density, -1.0 / 3.0);

        double rho13 = std::pow(density, 1.0 / 3.0);

        double expcr = std::exp(-c / rho13);

        double denom = 1.0 + d * rhom13;

        double omega = std::exp(-c * rhom13) / denom * std::pow(density, -11.0 / 3.0);

        double delta = rhom13 * (c + d / denom);

        double t2 = (47.0 - 7.0 * delta) * density_gradient / 18.0;

        double t5 = -2.0 / 3.0 * rho2 * density_gradient;

        double delta2 = -2.0 * pair_density;

        double gradR = 0.25 * density_gradient * (1.0 + delta2 / rho2);

        double t3 = -(2.5 - delta / 18.0) * (2.0 * gradR);

        double t1; double t4; double t6;
        // Real case
        if (pair_density <= 0)
        {
            double zeta  = std::sqrt(-2.0 * pair_density) / density;

            double gradI = 0.5 * density_gradient * zeta;

            double rhoa = 0.5 * density * (1.0 + zeta);

            double rhob = 0.5 * density * (1.0 - zeta);

            double grada2 = gradR + gradI;

            double gradb2 = gradR - gradI;

            t1 = std::pow(2.0, 11.0 / 3.0) * cf * (std::pow(rhoa, f83) + std::pow(rhob, f83));

            t4 = (11.0 - delta) / 9.0 * (rhoa * grada2 + rhob * gradb2) / density;

            t6 = ((2.0 / 3.0 * rho2 - rhoa * rhoa) * gradb2 + (2.0 / 3.0 * rho2 - rhob * rhob) * grada2);
        }
        // Imaginary case
        else
        {
            double eta    = std::sqrt(2.0 * pair_density) / density;

            double gradI = 0.5 * density_gradient * eta;

            // rho_a/b
            double R = (0.5*density) * std::sqrt(1 + eta*eta);

            double arg = std::atan(eta);

            // grad**2
            double Rgrad2 = std::sqrt(gradR*gradR + gradI*gradI);

            double arggrad2 = std::atan2(gradI,gradR);

            t1 = std::pow(2.0, 14.0 / 3.0) * cf * std::pow(R, f83) * std::cos(f83*arg);

            t4 = 2.0 * (11.0 - delta) / 9.0 * Rgrad2 * R * std::cos(arg + arggrad2)/density;

            t6 = f43 * rho2 * gradR - 2 * R*R * Rgrad2 * std::cos(-2.0*arg + arggrad2);
        }

        double rarb = 0.25 * density*density * (1.0 - delta2 / rho2);

        exc[g] = -a * (4.0 * rarb / (denom * density) + b * omega * (rarb * (t1 + t2 + t3 + t4) + t5 + t6)) / density;

        vrho[2 * g + 0] = 0.0; 
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_c
