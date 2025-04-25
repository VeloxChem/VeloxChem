//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "PairDensitySlater.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "MathConst.hpp"

namespace pdftslater {  // pdftslater namespace

void
compute_exc_vxc(const int np, const double* rho, double* exc, double* vrho)
{
    double onethird = 1.0 / 3.0;

    double twothird = 2.0 / 3.0;

    double fourthird = 4.0 / 3.0;

    double onesix = 1.0 / 6.0;

    double eightnine = 8.0 / 9.0;

    double eightthird = 8.0 / 3.0;

    double frg = -std::pow(6.0 / M_PI, onethird);

    double fre = frg * 0.75 / std::pow(2.0, fourthird);

    for (int g = 0; g < np; g++)
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

            double dxa = std::pow(1.0 + zeta, onethird);

            f_zeta = fxa;

            double fl_zeta = dxa;

            if (1.0 - zeta > 1.0e-16)
            {
                double fxb = std::pow(1.0 - zeta, fourthird);

                f_zeta += fxb;

                double dxb = std::pow(1.0 - zeta, onethird);

                fl_zeta -= dxb;
            }

            // dExc/d(rho)

            dexc_rho = fre * fourthird * rhothird * (f_zeta - zeta * fl_zeta);

            // dExc/d(Pi)

            if (pair_density < -1.0e-16)
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

            double pref = eightthird * std::pow(rs, -onethird) * eta / density;

            double df_zeta_drho = pref * (-eta * std::cos(theta) + std::sin(theta));

            // dExc/d(rho)

            dexc_rho = fre * rhothird * (fourthird * f_zeta + density * df_zeta_drho);

            // dExc/d(Pi)

            if (pair_density > 1.0e-16)
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
