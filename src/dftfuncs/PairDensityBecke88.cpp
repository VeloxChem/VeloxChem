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

#include "PairDensityBecke88.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftb88 {  // pdftb88 namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double fb = 0.0042; 

    double f43 = 4.0 / 3.0;

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

        double sig = sigma[3 * g + 0];

        // Li Manni's 2014 translation of gradients
        double delta2 = -2.0 * pair_density;

        double rho2 = std::pow(density, 2.0);

        double gradR = 0.25 * sig * (1.0 + delta2 / rho2);

        double fzeta = 0.0;
        // Real case
        if (pair_density <= 0)
        {
            double delta = std::sqrt(-2.0 * pair_density);

            double zeta = delta / density;

            double gradI = 0.5 * sig * zeta;

            double ra = std::pow(0.5*(density + delta), f43);

            double rb = std::pow(0.5*(density - delta), f43);

            double grada2 = gradR + gradI;

            double gradb2 = gradR - gradI;

            double grada = std::sqrt(std::abs(grada2)); // x asinh(x) is even function

            double gradb = std::sqrt(std::abs(gradb2));

            double xa = grada /ra;

            double xb = gradb /rb;

            double xa2 = grada2 / (ra*ra);

            double xb2 = gradb2 / (rb*rb);

            double fda = 1.0 + 6.0 * fb * xa * std::asinh(xa);

            double fdb = 1.0 + 6.0 * fb * xb * std::asinh(xb);

            double fea = ra * xa2 / fda;

            double feb = rb * xb2 / fdb;

            fzeta = fea + feb;
        }
        // Imaginary case
        else
        {
            double delta = std::sqrt(2.0 * pair_density);

            double eta = delta / density;

            double gradI = 0.5 * sig * eta;

            // rho**4/3
            double R = std::pow(0.5*density,f43) * std::pow(1.0 + eta*eta,0.5*f43);

            double arg = f43*std::atan(eta);

            // grad**2
            double Rgrad2 = std::sqrt(gradR*gradR + gradI*gradI);

            double arggrad2 = std::atan2(gradI,gradR);

            // x**2 = grad**2 / rho**8/3
            double Rx2 = Rgrad2/(R*R);

            double argx2 = arggrad2 - 2 * arg;

            double x2_R = Rx2 * std::cos(argx2);

            double x2_I = Rx2 * std::sin(argx2);

            // x = sqrt(x**2)
            double R_x = std::sqrt(Rx2);

            double x_R = R_x * std::cos(0.5*argx2);

            double x_I = R_x * std::sin(0.5*argx2);

            // sqrt(1 + x**2)
            double Rx2_1 = std::pow((1.0+x2_R)*(1.0+x2_R)+x2_I*x2_I, 0.25);

            double argx2_1 = 0.5*std::atan2(x2_I ,1+x2_R);

            // x + sqrt(1 + x**2)
            double ln_R = Rx2_1 * std::cos(argx2_1) + x_R;

            double ln_I = Rx2_1 * std::sin(argx2_1) + x_I;

            // arcsinh(x) = log(x + sqrt(1 + x**2))
            double arcsinh_R = 0.5 * std::log(ln_R*ln_R + ln_I*ln_I);

            double arcsinh_I = std::atan2(ln_I,ln_R);

            // fd = 1 + 6 * fb * x * arcsinh(x)
            double fd_R =  1.0 + 6.0 * fb * (arcsinh_R*x_R - arcsinh_I*x_I);

            double fd_I = 6.0 * fb * (arcsinh_R * x_I + arcsinh_I * x_R);

            double R_fd = std::sqrt(fd_R*fd_R + fd_I*fd_I);

            double arg_fd = std::atan2(fd_I, fd_R);

            // fzeta = rho**4/3 * x**2 / fd
            fzeta = 2.0 * R * Rx2 / R_fd * std::cos(arg + argx2 - arg_fd);
        }

        exc[g] = -fb * fzeta/density;

        vrho[2 * g + 0] = 0.0; 
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftb88
