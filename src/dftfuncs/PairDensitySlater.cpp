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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftslater {  // pdftslater namespace

void
compute_exc_vxc(const int32_t np, const double* rho, double* exc, double* vrho)
{
    double onethird = 1.0 / 3.0;

    double twothird = 2.0 / 3.0;

    double fourthird = 4.0 / 3.0;

    double onesix = 1.0 / 6.0;

    double eightnine = 8.0 / 9.0;

    double eightthird = 8.0 / 3.0;

    double frg = -std::pow(6.0 / mathconst::getPiValue(), onethird);

    double fre = frg * 0.75 / std::pow(2.0, fourthird);

    for (int32_t g = 0; g < np; g++)
    {
        double density = rho[2 * g + 0];

        if (density < 1.0e-16)
        {
            exc[g] = 0.0;

            vrho[2 * g + 0] = 0.0;

            vrho[2 * g + 1] = 0.0;

            continue;
        }

        double rhothird = std::pow(density, onethird);

        // double rhofourthird = std::pow(density, fourthird);

        double rhotwothird = std::pow(density, twothird);

        double pair_density = rho[2 * g + 1];

        double f_zeta = 0.0;

        double dexc_rho = 0.0;

        double dexc_pi = 0.0;

        // Real case
        if (pair_density <= 0)
        {
            // Energy

            double delta = sqrt(-2.0 * pair_density);

            double zeta = delta / density;

            double fxa = std::pow(1.0 + zeta, fourthird);

            double fxb = std::pow(1.0 - zeta, fourthird);

            f_zeta = fxa + fxb;

            // Derivatives

            double dxa = std::pow(1.0 + zeta, onethird);

            double dxb = std::pow(1.0 - zeta, onethird);

            double fl_zeta = dxa - dxb;

            // dExc/d(rho)

            dexc_rho = fre * fourthird * rhothird * (f_zeta - zeta * fl_zeta);

            // dExc/d(Pi)

            if (pair_density < -1.0e-12)
            {
                dexc_pi = -fre * fourthird * rhothird / delta * fl_zeta;
            }
            else
            {
                dexc_pi = -fre * eightnine / rhotwothird;
            }
        }
        // Imaginary case
        else
        {
            // Energy

            double delta = sqrt(2.0 * pair_density);

            double eta = delta / density;

            double rs = 1.0 + std::pow(eta, 2);

            double theta = fourthird * std::atan(eta);

            f_zeta = 2.0 * std::pow(rs, twothird) * std::cos(theta);

            // Derivatives

            double gr = std::pow(rs, onesix);

            double gtheta = onethird * std::atan(eta);

            double flim_eta = std::sin(gtheta);

            // dExc/d(rho)

            dexc_rho = fre * fourthird * rhothird * (f_zeta - eta * 2.0 * gr * flim_eta);

            // dExc/d(Pi)

            if (pair_density > 1.0e-12)
            {
                dexc_pi = -fre * eightthird * rhothird / delta * gr * flim_eta;
            }
            else
            {
                dexc_pi = -fre * eightnine / rhotwothird;
            }
        }

        exc[g] = fre * rhothird * f_zeta;

        // dE/drho
        vrho[2 * g + 0] = dexc_rho;
        // dE/dPi
        vrho[2 * g + 1] = dexc_pi;
    }
}

}  // namespace pdftslater
