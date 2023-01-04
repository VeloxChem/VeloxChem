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

#include "PairDensityPBE_C.hpp"

#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftpbe_c {  // pdftpbe_c namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double onethird  = 1.0 / 3.0;
    double twothird  = 2.0 / 3.0;
    double fourthird = 4.0 / 3.0;

    double frg = std::pow(0.75 / mathconst::getPiValue(), onethird);

    double two13 = std::pow(2.0, onethird);
    double two23 = two13 * two13;

    double three13 = std::pow(3.0, onethird);

    double pi13 = std::pow(mathconst::getPiValue(), onethird);

    double l       = 0.0716;
    double lambda  = 0.004235 * frg * two23 * 16.0;
    double lambda2 = lambda * lambda;

    double d2fact = 1.0 / 16.0 * pi13 / three13;

    double c = 1.709921;

    double t1 = 0.031091;
    double u1 = 0.21370;
    double v1 = 7.5957;
    double w1 = 3.5876;
    double x1 = 1.6382;
    double y1 = 0.49294;

    double t2 = 0.015545;
    double u2 = 0.20548;
    double v2 = 14.1189;
    double w2 = 6.1977;
    double x2 = 3.3662;
    double y2 = 0.62517;

    double t3 = 0.016887;
    double u3 = 0.11125;
    double v3 = 10.357;
    double w3 = 3.6231;
    double x3 = 0.88026;
    double y3 = 0.49671;

    double omega_fact = 1.0 / (2.0 * two13 - 2.0);

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

        double zeta4 = 4.0 * std::pow(pair_density, 2.0) / std::pow(density, 4.0);

        // Real case
        double omega;
        double u_ab;
        if (pair_density <= 0)
        {
            double delta = std::sqrt(-2.0 * pair_density);
            double zeta  = delta / density;

            omega = omega_fact * (std::pow(1.0 + zeta, fourthird) + std::pow(std::max(1.0 - zeta, 0.0), fourthird) - 2.0);

            u_ab = 0.5 * std::pow(1.0 + zeta, twothird) + 0.5 * std::pow(std::max(1.0 - zeta, 0.0), twothird);
        }
        // Imaginary case
        else
        {
            double delta  = std::sqrt(2.0 * pair_density);
            double eta    = delta / density;
            double r      = 1.0 + std::pow(eta, 2);
            double theta  = fourthird * std::atan(eta);
            double theta2 = twothird * std::atan(eta);

            omega = omega_fact * (2.0 * std::pow(r, twothird) * std::cos(theta) - 2.0);

            u_ab = std::pow(r, onethird) * std::cos(theta2);
        }

        double u_ab2 = u_ab * u_ab;

        double u_ab3 = u_ab2 * u_ab;

        // Compute epsilon
        double r = frg / std::pow(density, onethird);

        double r2 = r * r;

        double r12 = std::sqrt(r);

        double r32 = r * r12;

        double e1 = -2.0 * t1 * (u1 * r + 1.0) * log(0.5 / (t1 * (v1 * r12 + w1 * r + x1 * r32 + y1 * r2)) + 1.0);

        double e2 = -2.0 * t2 * (u2 * r + 1.0) * log(0.5 / (t2 * (v2 * r12 + w2 * r + x2 * r32 + y2 * r2)) + 1.0);

        double e3 = -2.0 * t3 * (u3 * r + 1.0) * log(0.5 / (t3 * (v3 * r12 + w3 * r + x3 * r32 + y3 * r2)) + 1.0);

        double epsilon = e1 - e3 * omega * (1.0 - zeta4) / c + (e2 - e1) * omega * zeta4;

        // Gradient-dependent terms
        double density_gradient = sigma[3 * g + 0];

        double d2 = d2fact * density_gradient / (u_ab2 * std::pow(density, 7.0 / 3.0));

        double d4 = d2 * d2;

        double N_ab = 2.0 * l / lambda / (exp(-2.0 * l / lambda2 * epsilon / u_ab3) - 1.0);

        exc[g] = 0.5 * lambda2 / l * u_ab3 * log(1.0 + 2.0 * l / lambda * (N_ab * d4 + d2) / (N_ab * d2 + N_ab * N_ab * d4 + 1.0)) + epsilon;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_c
