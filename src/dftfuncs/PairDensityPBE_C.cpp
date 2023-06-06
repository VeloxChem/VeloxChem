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

#include <algorithm>
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

        double dzeta4_drho = - 4.0 * zeta4 * (1.0 / density);

        double dzeta4_dpi = 8.0 * pair_density / std::pow(density, 4.0);

        double omega;
        double u_ab;
        double dom_drho;
        double dom_dpi;
        double dmu_drho;
        double dmu_dpi;

        // Real case
        if (pair_density <= 0)
        {
            double delta = std::sqrt(-2.0 * pair_density);

            double zeta  = delta / density;

            omega = omega_fact * (std::pow(1.0 + zeta, fourthird) + std::pow(1.0 - zeta, fourthird) - 2.0);

            u_ab = 0.5 * std::pow(1.0 + zeta, twothird) + 0.5 * std::pow(1.0 - zeta, twothird);

            double dom_dzeta = omega_fact * fourthird * (std::pow(1.0+zeta, onethird)-std::pow(1.0-zeta,onethird));

            double dmu_dzeta = onethird * (std::pow(1.0+zeta,-onethird)-std::pow(1.0-zeta,-onethird));

            double dzeta_drho = - zeta/density;

            dom_drho = dom_dzeta * dzeta_drho;

            dmu_drho = dmu_dzeta * dzeta_drho;

            if (pair_density < -1.0e-12)
            {
                double dzeta_dpi = - 1.0 /(delta * density);

                dom_dpi = dom_dzeta * dzeta_dpi;

                dmu_dpi= dmu_dzeta * dzeta_dpi;
            }
            else
            {
                dom_dpi = - omega_fact * fourthird * twothird / (density * density);

                dmu_dpi = onethird * twothird / (density * density);
            }

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

            double dom_dr = omega_fact * (4 / 3) * std::pow (r,(2/3)-1)*std::cos(theta);

            double dr_deta = 2 * eta;

            double deta_drho = - std::pow(2 * pair_density, 0.5) / std::pow(density, 2);

            double dmu_dr = (1/3)*std::pow(r,-(2/3))*std::cos(theta2);

            double dtheta2_deta = (2/3)*std::pow((std::pow(eta,2)+1),-1);

            double dom_dtheta = -2 * omega_fact * std::pow(r,(2/3)) * std::sin(theta);

            double dtheta_deta = 2 * dtheta2_deta;

            double dmu_dtheta2 = -1 * std::pow(r,(1/3))*std::sin(theta2);

            dom_drho = dom_dr * dr_deta * deta_drho + dom_dtheta * dtheta_deta * deta_drho;

            dmu_drho = dmu_dr * dr_deta * deta_drho + dmu_dtheta2 * dtheta2_deta * deta_drho;

            if (pair_density > 1.0e-12)
            {
                double deta_dpi = 1/(std::pow(2*pair_density,0.5)*density);

                dom_dpi = dom_dr * dr_deta * deta_dpi + dom_dtheta * dtheta_deta * deta_dpi;

                dmu_dpi = dmu_dr * dr_deta * deta_dpi + dmu_dtheta2 * dtheta2_deta * deta_dpi;
            }
            else
            {
                dom_dpi = - omega_fact * fourthird * twothird / (density * density);

                dmu_dpi = onethird * twothird / (density * density);
            }
        }

        double u_ab2 = u_ab * u_ab;

        double u_ab3 = u_ab2 * u_ab;

        // Compute epsilon
        double r = frg / std::pow(density, onethird);

        double r2 = r * r;

        double r12 = std::sqrt(r);

        double r32 = r * r12;

        double r3 = r * r * r;

        double dr_drho = -(1.0 / 3.0) * frg * std::pow(density,-(4.0 / 3.0));

        double h1 = -2.0 * t1 * (u1 * r + 1.0);

        double q1 = t1 * (v1 * r12 + w1 * r + x1 * r32 + y1 * r2);

        double g1 = std::log((0.5/q1) + 1.0);

        double e1 = h1 * g1;

        double dh1_dr = - 2.0 * t1 * u1;

        double dg1_dq1 = -0.5 / (std::pow(q1,2)+0.5*q1);

        double dq1_dr = t1*(0.5*v1*std::pow(r,-(1.0/2.0))+w1+x1*(3.0/2.0)*std::pow(r,0.5)+2.0*y1*r);

        double h2 = -2.0 * t2 * (u2 * r + 1.0);

        double q2 = t2 * (v2 * r12 + w2 * r + x2 * r32 + y2 * r2);

        double g2 = std::log((0.5 / q2) + 1.0);

        double e2 = h2 * g2;

        double dh2_dr = - 2.0 * t2 * u2;

        double dg2_dq2 = -0.5 / (std::pow(q2, 2) + 0.5 * q2);

        double dq2_dr = t2 * (0.5 * v2 * std::pow(r, -(1.0 / 2.0)) + w2 + x2 * (3.0 / 2.0) * std::pow(r,0.5) + 2.0 * y2 * r);

        double h3 = -2.0 * t3 * (u3 * r + 1.0);

        double q3 = t3 * (v3 * r12 + w3 * r + x3 * r32 + y3 * r2);

        double g3 = std::log((0.5 / q3) + 1.0);

        double e3 = h3 * g3;

        double dh3_dr = - 2.0 * t3 * u3;

        double dg3_dq3 = -0.5 / (std::pow(q3, 2) + 0.5 * q3);

        double dq3_dr = t3 * (0.5 * v3 * std::pow(r, -(1.0 / 2.0)) + w3 + x3 * (3.0 / 2.0) * std::pow(r,0.5) + 2.0 * y3 * r);

        double epsilon = e1 - e3 * omega * (1.0 - zeta4) / c + (e2 - e1) * omega * zeta4;

        double der1_drho = dr_drho*(dh1_dr*g1 + dg1_dq1*dq1_dr*h1);

        double der2_drho = dr_drho*(dh2_dr*g2 + dg2_dq2*dq2_dr*h2);

        double der3_drho = dr_drho*(dh3_dr*g3 + dg3_dq3*dq3_dr*h3);

        double deps_de1 = 1-omega*zeta4;

        double deps_de2 = omega*zeta4;

        double deps_de3 = -omega*(1-zeta4);

        double deps_dom = zeta4*(e2-e1)-(e3/c)*(1+zeta4);

        double deps_dzeta4 = omega*(e2-e1-(e3/c));

        double depsilon_drho = deps_de1 * der1_drho + deps_de2 * der2_drho + deps_de3 * der3_drho + deps_dom *  dom_drho + deps_dzeta4 * dzeta4_drho;

        double depsilon_dpi =  deps_dom * dom_dpi + deps_dzeta4 * dzeta4_dpi;

        // Gradient-dependent terms
        double density_gradient = sigma[3 * g + 0];

        double d2 = d2fact * density_gradient / (u_ab2 * std::pow(density, 7.0 / 3.0));

        double dd_dsigma = d2fact / (u_ab2 * std::pow(density, 7.0 / 3.0));

        double dd_drho = -d2 * (2*(1.0/u_ab)*dmu_drho +(7.0/3.0)*(1/density));

        double dd_dpi = -d2 * (2/u_ab) * dmu_dpi;

        double d4 = d2 * d2;

        double G = -2.0 * l / lambda2 * epsilon / u_ab3;

        double hx = exp(G);

        double N_ab =  (2.0 * l / lambda)* 1 / ( hx - 1.0);

        double N2 = N_ab * N_ab;

        double f2 =-2.0 * l / lambda2;

        double f = 2.0 * l / lambda;

        double Q = 0.5 * lambda2/l *u_ab3;

        double P = std::log(1.0 + f * (N_ab * d4 + d2) / (N_ab * d2 + N2 * d4 + 1.0));

        exc[g] = Q * P + epsilon;

        double dQ_dmu = 0.5 * 3.0 *(lambda2/l)*u_ab2;

        double dQ_drho = dQ_dmu * dmu_drho;

        double dQ_dpi = dQ_dmu * dmu_dpi;

        double dG_drho = (2*l/lambda2)*((3*epsilon/(u_ab3*u_ab))*dmu_drho - depsilon_drho / u_ab3);

        double dG_dpi = (2*l/lambda2)*((3*epsilon/(u_ab3*u_ab))*dmu_dpi - depsilon_dpi / u_ab3);

        double U = exp(G)-1;

        double M = 1/U;

        double dM_dU = -1/std::pow(U,2);

        double dP_dN = -f*(d4*d2*N_ab*(d2*N_ab+2))/((N_ab * d2 + N2 * d4 + 1.0)*(f* (N_ab * d4 + d2) + N2 * d4 + d2*N_ab+1));

        double dP_dd =  1/((std::pow(d2*N_ab,2)+ d2*N_ab + 1)*(f*d2*(d2*N_ab+1) + std::pow(d2*N_ab,2) + d2*N_ab + 1));

        dP_dd = dP_dd* f*(2*d2*N_ab + 1);

        double dN_drho = f * dM_dU * dG_drho * exp(G);

        double dN_dpi = f * dM_dU * dG_dpi * exp(G);

        double dP_drho = dP_dN * dN_drho + dP_dd * dd_drho;

        double dP_dpi = dP_dN * dN_dpi + dP_dd * dd_dpi;

        vrho[2 * g + 0] = (P * dQ_drho + Q * dP_drho + depsilon_drho)*density + exc[g];
        vrho[2 * g + 1] = (P * dQ_dpi + Q * dP_dpi + depsilon_dpi)*density;

        vsigma[3 * g + 0] = density*(Q * dP_dd * dd_dsigma);
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }
}

}  // namespace pdftpbe_c
