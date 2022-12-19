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

#include "PairDensitySlater.hpp"

#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftslater {  // pdftslater namespace

void
compute_exc_vxc(const int32_t np, const double* rho, double* exc, double* vrho)
{
    double frg = -std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0);

    double fre = 0.75 * frg;

    double fp = 1.0 / 3.0;

    double twothird = 2.0 / 3.0;

    double fourthird = 4.0 / 3.0;

    fre /= std::pow(2.0, fourthird);

    for (int32_t g = 0; g < np; g++)
    {
        double density = rho[2 * g + 0];

        if (density < 1.0e-8)
        {
            exc[g] = 0.0;

            vrho[2 * g + 0] = 0.0;

            vrho[2 * g + 1] = 0.0;

            continue;
        }

        double pair_density = rho[2 * g + 1];

        double f_zeta;

        // Real case
        if (pair_density <= 0)
        {
            double delta = sqrt(-2.0 * pair_density);

            double zeta = delta / density;

            double fxa = std::pow(1.0 + zeta, fourthird);

            double fxb = std::pow(std::max(1.0 - zeta, 0.0), fourthird);

            f_zeta = fxa + fxb;
        }
        // Imaginary case
        else
        {
            double delta = sqrt(2.0 * pair_density);

            double zeta = delta / density;

            double r = sqrt(1.0 + std::pow(zeta, 2));

            double theta = 4.0 / 3.0 * std::atan(zeta);

            f_zeta = 2.0 * std::pow(r, 4.0 / 3.0) * std::cos(theta);
        }

        double rhothird = std::pow(density, fp);

        exc[g] = fre * rhothird * f_zeta;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;
    }
}

}  // namespace pdftslater
