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
    double frg = - std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0);

    double fre = 0.75 * frg;

    double fp = 1.0 / 3.0;

    double twothird = 2.0 / 3.0;

    double fourthird = 4.0 / 3.0;

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

        if (pair_density <= 0)
        {
            double delta = sqrt(-2.0 * pair_density);

            double rhoa = 0.5 * (density + delta);

            double rhob = 0.5 * (density - delta);

            double fxa = std::pow(rhoa, fp);

            double fxb = std::pow(rhob, fp);

            exc[g] = fre * (rhoa * fxa +  rhob * fxb);
        }
        else
        {
            double delta = sqrt(2.0 * pair_density);

            double r2 = 0.25 * std::pow(density, 2) + 0.5 * pair_density;

            double theta = fourthird * std::atan(delta/density);

            exc[g] = 2.0 * fre * std::pow(r2, twothird) * std::cos(theta);
        }
        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;
    }
}

}  // namespace pdftslater
