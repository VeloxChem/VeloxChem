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

#include "MathConst.hpp"

namespace pdftpbe_c {  // pdftpbe_c namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
    double onethird = 1.0 / 3.0;
    double fourthird = 4.0 / 3.0;

    double frg = std::pow(0.75 / mathconst::getPiValue(), onethird);
    double frg2 = frg * frg; 
    double frg3 = sqrt(frg);

    double t2 = std::pow(3.0, onethird);
    double t4 = std::pow(mathconst::getPiValue(), onethird);
    double pi_43 = std::pow(mathconst::getPiValue(),fourthird);
    double t5 = 1.0 / t4;
    double t16 = std::pow(2.0,onethird);
    double t36 = 1.0/t2;

    double t8 = t2 * t2;
    double t10 = t4 * t4;
    double t11 = t5 * t5;
    double t31 = t16 * t16;
    double t41 = t36 * t36;

    double t13 = 1.732050807568877 * 0.5641895835477563;

    double t23 = 1.0/(2.0*t16-2.0);

    for (int32_t g = 0; g < np; g++)
    {
        double density = rho[2 * g + 0];  //t1

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

        double density2 = density * density;
        double inv_density = 1.0/density;
        double inv_density2 = inv_density * inv_density;
        double t15 = 1.0/sqrt(density);

        double pair_density = rho[2 * g + 1];
        double zeta4 = 4.0 * std::pow(pair_density,2.0)/std::pow(density,4.0);

        double t6 = 1.0/ std::pow(density,onethird);
        double t12 = t6 * t6;
        double t20 = sqrt(t6);

        double t7 = 0.2137*frg*t6+1.0;
        double t21 = log(16.0818243221511/(7.5957*frg3*t20+3.5876*
        frg*t6+0.8191*t13*t15+0.49294*frg2*t12)+1.0);
        double t22 = -0.062182*t7*t21;
        double reduced_sigma = sigma[3 * g + 0] * inv_density;
        double reduced_sigma2 = reduced_sigma * reduced_sigma;
        // Real case
        double t28;
        double t32;
        if (pair_density <= 0)
        {
            double delta = sqrt(-2.0 * pair_density);
            double zeta = delta * inv_density;
            t28 = std::pow(1.0+zeta,fourthird)+std::pow(1.0-zeta,fourthird)-2.0;
            double rhoa = 0.5* (density + delta);
            double rhob = 0.5* (density - delta);
            t32 = t31* std::pow(rhob,.6666666666666666)+t31*std::pow(rhoa,.6666666666666666);
        }
        // Imaginary case
        else
        {
            t28 = 0.0;
            t32 = t31 * t16 * std::pow(density, .6666666666666666);
        }
        double t29 = 0.0197517897025652*log(29.60857464321668/(10.357*
        frg3*t20+3.6231*frg*t6+0.44013*t13*t15+0.49671*
        frg2*t12)+1.0)*t23*(1.0-zeta4)*t28*(0.11125*frg*t6+1.0);

        double t30 = t23*zeta4*t28*(0.062182*t7*t21-0.03109*log(32.1646831778707/
        (14.1189*frg3*t20+6.1977*frg*t6+1.6831*t13*
        t15+0.62517*frg2*t12)+1.0)*(0.20548*frg*t6+1.0));

        double t33 = std::pow(t32, 3.0);
        double t38 = 1.0/std::pow(t32,2.0);
        double t40 = t38 * t38;
        double t42 = 1.0/exp(249.5089969598932*t10*(t30+t29+t22)*density2*t41/t33) -1.0;
        double t43 = 1.0/t42;
        double t42_2 = t42 * t42;

        exc[g] = 0.0040078715083799*t11*t33*inv_density2*log(2.113341204250295*
        t36*t4*(0.1383178583708611*reduced_sigma2*t40*t43+0.25*t36*t4*reduced_sigma*
        t38)/(0.5283353010625738*t41*t10*reduced_sigma*t38*t43+0.0930460634496268*
        pi_43*t36*reduced_sigma2*t40/t42_2+1.0)+1.0)*t8+t30+t29+t22;


        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_c
