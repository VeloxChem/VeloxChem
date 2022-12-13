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

#include "PairDensityVWN.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace pdftvwn {  // pdftvwn namespace

void
compute_exc_vxc(const int32_t np, const double* rho, double* exc, double* vrho)
{
    // paramagnetic fitting factors

    const double spinpolf = 1.92366105093154;

    const double fourthree = 1.333333333333333;

    double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;

    double pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578, px0_f = -0.7432940;

    // Paramagnetic

    double Q = std::sqrt(4.0 * pc - pb * pb);

    double XF0I = px0 * px0 + pb * px0 + pc;

    double YF0I = Q / (pb + 2.0 * px0);

    double B = px0 / XF0I;

    double C = XF0I * YF0I;

    double ACON = B * pb - 1.0;

    double BCON = 2.0 * ACON + 2.0;

    double CCON = 2.0 * pb * (1.0 / Q - px0 / C);

    // Ferromagnetic

    double Q_f = std::sqrt(4.0 * pc_f - pb_f * pb_f);

    double XF0I_f = px0_f * px0_f + pb_f * px0_f + pc_f;

    double YF0I_f = Q_f / (pb_f + 2.0 * px0_f);

    double B_f = px0_f / XF0I_f;

    double C_f = XF0I_f * YF0I_f;

    double ACON_f = B_f * pb_f - 1.0;

    double BCON_f = 2.0 * ACON_f + 2.0;

    double CCON_f = 2.0 * pb_f * (1.0 / Q_f - px0_f / C_f);

    // various prefactor

    double f16 = -1.0 / 6.0;

    double f76 = -7.0 / 6.0;

    double DCRS = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);

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

        double x = DCRS * std::pow(density, f16);

        double xrho = DCRS * f16 * std::pow(density, f76);

        double f_zeta;

        double f_zet1;

        // Real case
        if (pair_density <= 0)
        {
            double delta = sqrt(-2.0 * pair_density);

            double zeta = delta / density;

            f_zeta = spinpolf * (std::pow(1.0 + zeta, fourthree) + std::pow(1.0 - zeta, fourthree) - 2.0);

            f_zet1 = spinpolf * 4.0 / 3.0 * (std::pow(1.0 + zeta, 1.0 / 3.0) - std::pow(1.0 - zeta, 1.0 / 3.0));
        }
        // Imaginary case
        else
        {
            double delta = sqrt(2.0 * pair_density);

            double zeta = delta / density;

            double r = sqrt(1.0 + std::pow(zeta, 2));

            double theta = 4.0 / 3.0 * std::atan(zeta);

            f_zeta = spinpolf * (2.0 * std::pow(r, 4.0 / 3.0) * std::cos(theta) - 2.0);

            f_zet1 = 0.0;
        }

        double xf = x * x + pb * x + pc;

        double xf_f = x * x + pb_f * x + pc_f;

        double xfx = 2.0 * x + pb;

        double xfx_f = 2.0 * x + pb_f;

        double yf = Q / xfx;

        double yf_f = Q_f / xfx_f;

        double e1 = 2.0 * std::log(x) + ACON * std::log(xf) - BCON * std::log(x - px0) + CCON * std::atan(yf);

        double e1_f = 2.0 * std::log(x) + ACON_f * std::log(xf_f) - BCON_f * std::log(x - px0_f) + CCON_f * std::atan(yf_f);

        double ep_p0 = 0.5 * pa * e1;

        double ep_f0 = 0.5 * pa_f * e1_f;

        // ep_p1

        double ex1 = 2.0 / x + ACON * xfx / xf - BCON / (x - px0) - CCON * (2.0 * yf / xfx) / (1.0 + yf * yf);

        double ex1_f = 2.0 / x + ACON_f * xfx_f / xf_f - BCON_f / (x - px0_f) - CCON_f * (2.0 * yf_f / xfx_f) / (1.0 + yf_f * yf_f);

        double ep_p1 = 0.5 * pa * (e1 + density * ex1 * xrho);

        double ep_f1 = 0.5 * pa_f * (e1_f + density * ex1_f * xrho);

        // Potential

        double ef0 = ep_f0 - ep_p0;

        double ef1 = ep_f1 - ep_p1;

        double vcfp = f_zeta * ef1;

        double delta = f_zet1 * ef0;

        double delta_en = f_zeta * ef0;

        exc[g] = ep_p0 + delta_en;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;
    }
}

}  // namespace pdftvwn
