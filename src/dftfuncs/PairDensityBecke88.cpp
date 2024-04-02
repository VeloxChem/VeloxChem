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
        // No explicit dependence
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;

        double density = rho[2 * g + 0];

        if (density < 1.0e-16)
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
        double ddelta2_dpi = -2;

        double rho2 = std::pow(density, 2.0);
        double drho2_drho = 2*density;

        double zeta2 = delta2/rho2;
        double dzeta2_drho = -delta2*drho2_drho/pow(rho2, 2);
        double dzeta2_dpi = ddelta2_dpi/rho2;

        double gradR = 0.25 * sig * (1.0 + zeta2);
        double dgradR_dzeta2 = 0.25 *sig;
        double dgradR_drho = dgradR_dzeta2 * dzeta2_drho;
        double dgradR_dpi = dgradR_dzeta2*dzeta2_dpi;
        double dgradR_dsig = (1.0/4.0)*zeta2 + 1.0/4.0;

        double fzeta = 0.0;
        double dfzeta_drho = 0.0;
        double dfzeta_dpi = 0.0;
        double dfzeta_dsig = 0.0;

        // Real case
        if (pair_density < -1.0e-30)
        {
            double zeta = std::sqrt(-2.0 * pair_density) / density;
            double dzeta_drho = - sqrt(-2.0 * pair_density)/pow(density, 2);
            double dzeta_dpi = (1.0/2.0)* sqrt(-2.0 * pair_density)/(density*pair_density);

            double gradI = 0.5 * sig * zeta;
            double dgradI_dzeta = (1.0/2.0)*sig;
            double dgradI_drho = dgradI_dzeta*dzeta_drho;
            double dgradI_dpi = dgradI_dzeta*dzeta_dpi;
            double dgradI_dsig = (1.0/2.0)*zeta;

            double ra = std::pow(0.5*density * (1+ zeta), f43);
            double dra_dzeta = f43*pow((1.0/2.0)*density*(zeta + 1), f43)/(zeta + 1);
            double dra_drho = dra_dzeta*dzeta_drho + 2*f43*pow((1.0/2.0)*density*(zeta + 1), f43)*((1.0/2.0)*zeta + 1.0/2.0)/(density*(zeta + 1));
            double dra_dpi = dra_dzeta*dzeta_dpi;

            double rb = std::pow(0.5*density* (1 - zeta), f43);
            double drb_dzeta = 0.0;
            double drb_drho = 0.0;
            if (rb > 1.0e-16)
            {
                drb_dzeta = -f43*pow((1.0/2.0)*density*(1 - zeta), f43)/(1 - zeta);
                drb_drho = drb_dzeta*dzeta_drho + 2*f43*pow((1.0/2.0)*density*(1 - zeta), f43)*(1.0/2.0 - 1.0/2.0*zeta)/(density*(1 - zeta));
            }
            double drb_dpi = drb_dzeta*dzeta_dpi;

            double grada2 = gradR + gradI;
            double dgrada2_drho = dgradR_drho + dgradI_drho;
            double dgrada2_dpi = dgradR_dpi + dgradI_dpi;
            double dgrada2_dsig = dgradR_dsig + dgradI_dsig;

            double gradb2 = std::max(0.0, gradR - gradI);
            double dgradb2_drho = dgradR_drho - dgradI_drho;
            double dgradb2_dpi = dgradR_dpi - dgradI_dpi;
            double dgradb2_dsig = dgradR_dsig - dgradI_dsig;

            double grada = std::sqrt(grada2);
            double dgrada_dgrada2 = (1.0/2.0)/sqrt(grada2);
            double dgrada_drho = dgrada2_drho*dgrada_dgrada2;
            double dgrada_dpi = dgrada2_dpi*dgrada_dgrada2;
            double dgrada_dsig = dgrada2_dsig*dgrada_dgrada2;

            double gradb = std::sqrt(gradb2);
            double dgradb_dgradb2 = 0.0;
            if (gradb2 > 1.0e-16)
            {
                dgradb_dgradb2 = (1.0/2.0)/sqrt(gradb2);
            }
            double dgradb_drho = dgradb2_drho*dgradb_dgradb2;
            double dgradb_dpi = dgradb2_dpi*dgradb_dgradb2;
            double dgradb_dsig = dgradb2_dsig*dgradb_dgradb2;

            double xa = grada /ra;
            double dxa_dra = -grada/pow(ra, 2);
            double dxa_dgrada = 1.0/ra;
            double dxa_drho = dgrada_drho*dxa_dgrada + dra_drho*dxa_dra;
            double dxa_dpi = dgrada_dpi*dxa_dgrada + dra_dpi*dxa_dra;
            double dxa_dsig = dgrada_dsig*dxa_dgrada;

            double xa2 = grada2 / (ra*ra);
            double dxa2_dra = -2*grada2/pow(ra, 3);
            double dxa2_dgrada2 = pow(ra, -2);
            double dxa2_drho = dgrada2_drho*dxa2_dgrada2 + dra_drho*dxa2_dra;
            double dxa2_dpi = dgrada2_dpi*dxa2_dgrada2 + dra_dpi*dxa2_dra;
            double dxa2_dsig = dgrada2_dsig*dxa2_dgrada2;

            double xb = 0.0;
            double dxb_drb = 0.0;
            double dxb_dgradb = 0.0;
            if (rb > 1.0e-16)
            {
                xb = gradb /rb;
                dxb_drb = -gradb/pow(rb, 2);
                dxb_dgradb = 1.0/rb;
            }
            double dxb_drho = dgradb_drho*dxb_dgradb + drb_drho*dxb_drb;
            double dxb_dpi = dgradb_dpi*dxb_dgradb + drb_dpi*dxb_drb;
            double dxb_dsig = dgradb_dsig*dxb_dgradb;

            double xb2 = 0.0;
            double dxb2_dgradb2 = 0.0;
            double dxb2_drb = 0.0;
            if (rb > 1.0e-16)
            {
                xb2 = gradb2 / (rb*rb);
                dxb2_dgradb2 = pow(rb, -2);
                dxb2_drb = -2*gradb2/pow(rb, 3);
            }
            double dxb2_drho = dgradb2_drho*dxb2_dgradb2 + drb_drho*dxb2_drb;
            double dxb2_dpi = dgradb2_dpi*dxb2_dgradb2 + drb_dpi*dxb2_drb;
            double dxb2_dsig = dgradb2_dsig*dxb2_dgradb2;

            double fda = 1.0 + 6.0 * fb * xa * std::asinh(xa);
            double dfda_dxa = 6*fb*xa/sqrt(pow(xa, 2) + 1) + 6*fb*asinh(xa);
            double dfda_drho = dfda_dxa*dxa_drho;
            double dfda_dpi = dfda_dxa*dxa_dpi;
            double dfda_dsig = dfda_dxa*dxa_dsig;

            double fdb = 1.0 + 6.0 * fb * xb * std::asinh(xb);
            double dfdb_dxb = 6*fb*xb/sqrt(pow(xb, 2) + 1) + 6*fb*asinh(xb);
            double dfdb_drho = dfdb_dxb*dxb_drho;
            double dfdb_dpi = dfdb_dxb*dxb_dpi;
            double dfdb_dsig = dfdb_dxb*dxb_dsig;

            double fea = ra * xa2 / fda;
            double dfea_dra = xa2/fda;
            double dfea_dxa2 = ra/fda;
            double dfea_dfda = -ra*xa2/pow(fda, 2);
            double dfea_drho = dfda_drho*dfea_dfda + dfea_dra*dra_drho + dfea_dxa2*dxa2_drho;
            double dfea_dpi = dfda_dpi*dfea_dfda + dfea_dra*dra_dpi + dfea_dxa2*dxa2_dpi;
            double dfea_dsig = dfda_dsig*dfea_dfda + dfea_dxa2*dxa2_dsig;

            double feb = rb * xb2 / fdb;
            double dfeb_dfdb = -rb*xb2/pow(fdb, 2);
            double dfeb_dxb2 = rb/fdb;
            double dfeb_drb = xb2/fdb;
            double dfeb_drho = dfdb_drho*dfeb_dfdb + dfeb_drb*drb_drho + dfeb_dxb2*dxb2_drho;
            double dfeb_dpi = dfdb_dpi*dfeb_dfdb + dfeb_drb*drb_dpi + dfeb_dxb2*dxb2_dpi;
            double dfeb_dsig = dfdb_dsig*dfeb_dfdb + dfeb_dxb2*dxb2_dsig;

            fzeta = fea + feb;
            dfzeta_drho = dfea_drho + dfeb_drho;
            dfzeta_dpi = dfea_dpi + dfeb_dpi;
            dfzeta_dsig = dfea_dsig + dfeb_dsig;
        }
        // pi=0 special case
        else if (pair_density < 1.0e-30)
        {
            double xa = 0.5 * sqrt(sig) * pow(0.5 * density, -f43);
            double dxa_drho = -f43*cbrt(2)*sqrt(sig)/pow(density, 7.0/3.0);
            double dxa_dsig = (1.0/2.0)*cbrt(2)/(pow(density, 4.0/3.0)*sqrt(sig));

            double fda = 1.0 + 6.0 * fb * xa * std::asinh(xa);
            double dfda_dxa = 6*fb*xa/sqrt(pow(xa, 2) + 1) + 6*fb*asinh(xa);
            double dfda_drho = dfda_dxa*dxa_drho;
            double dfda_dsig = dfda_dxa*dxa_dsig;

            fzeta = 0.5 * sig * pow(0.5 * density, -f43)  / fda;
            double dfzeta_dfda = -cbrt(2)*sig/(pow(density, 4.0/3.0)*pow(fda, 2));
            dfzeta_drho = dfda_drho*dfzeta_dfda - 4.0/3.0*cbrt(2)*sig/(pow(density, 7.0/3.0)*fda);
            dfzeta_dpi = 0;
            dfzeta_dsig = dfda_dsig*dfzeta_dfda + cbrt(2)/(pow(density, 4.0/3.0)*fda);
        }
        // Imaginary case
        else
        {
            double eta = std::sqrt(2.0 * pair_density) / density;
            double deta_drho = - eta / density;
            double deta_dpi = 1.0 / (density * sqrt(2*pair_density));

            double gradI = 0.5 * sig * eta;
            double dgradI_deta = 0.5 * sig;
            double dgradI_drho = deta_drho * dgradI_deta;
            double dgradI_dpi = deta_dpi * dgradI_deta;
            double dgradI_dsig = 0.5 * eta;

            // rho**4/3
            double R = std::pow(0.5*density,f43) * std::pow(1.0 + eta*eta,0.5*f43);
            double dR_deta = (1.0/3.0)*pow(2, 2.0/3.0)*pow(density, 4.0/3.0)*eta/cbrt(pow(eta, 2) + 1);
            double dR_drho = dR_deta*deta_drho + (1.0/3.0)*pow(2, 2.0/3.0)*cbrt(density)*pow(pow(eta, 2) + 1, 2.0/3.0);
            double dR_dpi = dR_deta*deta_dpi;

            double argr = f43*std::atan(eta);
            double dargr_deta = (4.0/3.0)/(pow(eta, 2) + 1);
            double dargr_drho = dargr_deta*deta_drho;
            double dargr_dpi = dargr_deta*deta_dpi;

            // grad**2
            double Rgrad2 = std::sqrt(gradR*gradR + gradI*gradI);
            double dRgrad2_dgradR = gradR/Rgrad2;
            double dRgrad2_dgradI = gradI/Rgrad2;
            double dRgrad2_drho = dRgrad2_dgradI*dgradI_drho + dRgrad2_dgradR*dgradR_drho;
            double dRgrad2_dpi = dRgrad2_dgradI*dgradI_dpi + dRgrad2_dgradR*dgradR_dpi;
            double dRgrad2_dsig = dRgrad2_dgradI*dgradI_dsig + dRgrad2_dgradR*dgradR_dsig;

            double arggrad2 = std::atan2(gradI,gradR);
            double darggrad2_dgradR = -gradI/(pow(gradI, 2) + pow(gradR, 2));
            double darggrad2_dgradI = gradR/(pow(gradI, 2) + pow(gradR, 2));
            double darggrad2_drho = darggrad2_dgradI*dgradI_drho + darggrad2_dgradR*dgradR_drho;
            double darggrad2_dpi = darggrad2_dgradI*dgradI_dpi + darggrad2_dgradR*dgradR_dpi;
            double darggrad2_dsig = darggrad2_dgradI*dgradI_dsig + darggrad2_dgradR*dgradR_dsig;

            // x**2 = grad**2 / rho**8/3
            double Rx2 = Rgrad2/(R*R);
            double dRx2_dRgrad2 = pow(R, -2);
            double dRx2_dR = -2*Rgrad2/pow(R, 3);
            double dRx2_drho = dR_drho*dRx2_dR + dRgrad2_drho*dRx2_dRgrad2;
            double dRx2_dpi = dR_dpi*dRx2_dR + dRgrad2_dpi*dRx2_dRgrad2;
            double dRx2_dsig = dRgrad2_dsig*dRx2_dRgrad2;

            double argx2 = arggrad2 - 2 * argr;
            double dargx2_drho = darggrad2_drho - 2.0 * dargr_drho;
            double dargx2_dpi = darggrad2_dpi - 2.0 * dargr_dpi;
            double dargx2_dsig = darggrad2_dsig;

            double x2_R = Rx2 * std::cos(argx2);
            double dx2_R_dargx2 = -Rx2*sin(argx2);
            double dx2_R_dRx2 = cos(argx2);
            double dx2_R_drho = dRx2_drho*dx2_R_dRx2 + dargx2_drho*dx2_R_dargx2;
            double dx2_R_dpi = dRx2_dpi*dx2_R_dRx2 + dargx2_dpi*dx2_R_dargx2;
            double dx2_R_dsig = dRx2_dsig*dx2_R_dRx2 + dargx2_dsig*dx2_R_dargx2;

            double x2_I = Rx2 * std::sin(argx2);
            double dx2_I_dargx2 = Rx2*cos(argx2);
            double dx2_I_dRx2 = sin(argx2);
            double dx2_I_drho = dRx2_drho*dx2_I_dRx2 + dargx2_drho*dx2_I_dargx2;
            double dx2_I_dpi = dRx2_dpi*dx2_I_dRx2 + dargx2_dpi*dx2_I_dargx2;
            double dx2_I_dsig = dRx2_dsig*dx2_I_dRx2 + dargx2_dsig*dx2_I_dargx2;

            // x = sqrt(x**2)
            double R_x = std::sqrt(Rx2);
            double dR_x_dRx2 = (1.0/2.0) / R_x;
            double dR_x_drho = dR_x_dRx2*dRx2_drho;
            double dR_x_dpi = dR_x_dRx2*dRx2_dpi;
            double dR_x_dsig = dR_x_dRx2*dRx2_dsig;

            double x_R = R_x * std::cos(0.5*argx2);
            double dx_R_dargx2 = -1.0/2.0*R_x*sin((1.0/2.0)*argx2);
            double dx_R_dR_x = cos((1.0/2.0)*argx2);
            double dx_R_drho = dR_x_drho*dx_R_dR_x + dargx2_drho*dx_R_dargx2;
            double dx_R_dpi = dR_x_dpi*dx_R_dR_x + dargx2_dpi*dx_R_dargx2;
            double dx_R_dsig = dR_x_dsig*dx_R_dR_x + dargx2_dsig*dx_R_dargx2;

            double x_I = R_x * std::sin(0.5*argx2);
            double dx_I_dargx2 = (1.0/2.0)*R_x*cos((1.0/2.0)*argx2);
            double dx_I_dR_x = sin((1.0/2.0)*argx2);
            double dx_I_drho = dR_x_drho*dx_I_dR_x + dargx2_drho*dx_I_dargx2;
            double dx_I_dpi = dR_x_dpi*dx_I_dR_x + dargx2_dpi*dx_I_dargx2;
            double dx_I_dsig = dR_x_dsig*dx_I_dR_x + dargx2_dsig*dx_I_dargx2;

            // sqrt(1 + x**2)
            double Rx2_1 = std::pow((1.0+x2_R)*(1.0+x2_R)+x2_I*x2_I, 0.25);
            double dRx2_1_dx2_R = ((1.0/2.0)*x2_R + 1.0/2.0)/pow(pow(x2_I, 2) + pow(x2_R + 1, 2), 3.0/4.0);
            double dRx2_1_dx2_I = (1.0/2.0)*x2_I/pow(pow(x2_I, 2) + pow(x2_R + 1, 2), 3.0/4.0);
            double dRx2_1_drho = dRx2_1_dx2_I*dx2_I_drho + dRx2_1_dx2_R*dx2_R_drho;
            double dRx2_1_dpi = dRx2_1_dx2_I*dx2_I_dpi + dRx2_1_dx2_R*dx2_R_dpi;
            double dRx2_1_dsig = dRx2_1_dx2_I*dx2_I_dsig + dRx2_1_dx2_R*dx2_R_dsig;

            double argx2_1 = 0.5*std::atan2(x2_I ,1+x2_R);
            double dargx2_1_dx2_R = -1.0/2.0*x2_I/(pow(x2_I, 2) + pow(x2_R + 1, 2));
            double dargx2_1_dx2_I = (1.0/2.0)*(x2_R + 1)/(pow(x2_I, 2) + pow(x2_R + 1, 2));
            double dargx2_1_drho = dargx2_1_dx2_I*dx2_I_drho + dargx2_1_dx2_R*dx2_R_drho;
            double dargx2_1_dpi = dargx2_1_dx2_I*dx2_I_dpi + dargx2_1_dx2_R*dx2_R_dpi;
            double dargx2_1_dsig = dargx2_1_dx2_I*dx2_I_dsig + dargx2_1_dx2_R*dx2_R_dsig;

            // x + sqrt(1 + x**2)
            double ln_R = Rx2_1 * std::cos(argx2_1) + x_R;
            double dln_R_dargx2_1 = -Rx2_1*sin(argx2_1);
            double dln_R_dRx2_1 = cos(argx2_1);
            double dln_R_drho = dRx2_1_drho*dln_R_dRx2_1 + dargx2_1_drho*dln_R_dargx2_1 + dx_R_drho;
            double dln_R_dpi = dRx2_1_dpi*dln_R_dRx2_1 + dargx2_1_dpi*dln_R_dargx2_1 + dx_R_dpi;
            double dln_R_dsig = dRx2_1_dsig*dln_R_dRx2_1 + dargx2_1_dsig*dln_R_dargx2_1 + dx_R_dsig;

            double ln_I = Rx2_1 * std::sin(argx2_1) + x_I;
            double dln_I_dargx2_1 = Rx2_1*cos(argx2_1);
            double dln_I_dRx2_1 = sin(argx2_1);
            double dln_I_drho = dRx2_1_drho*dln_I_dRx2_1 + dargx2_1_drho*dln_I_dargx2_1 + dx_I_drho;
            double dln_I_dpi = dRx2_1_dpi*dln_I_dRx2_1 + dargx2_1_dpi*dln_I_dargx2_1 + dx_I_dpi;
            double dln_I_dsig = dRx2_1_dsig*dln_I_dRx2_1 + dargx2_1_dsig*dln_I_dargx2_1 + dx_I_dsig;

            // arcsinh(x) = log(x + sqrt(1 + x**2))
            double arcsinh_R = 0.5 * std::log(ln_R*ln_R + ln_I*ln_I);
            double darcsinh_R_dln_R = ln_R/(pow(ln_I, 2) + pow(ln_R, 2));
            double darcsinh_R_dln_I = ln_I/(pow(ln_I, 2) + pow(ln_R, 2));
            double darcsinh_R_drho = darcsinh_R_dln_I*dln_I_drho + darcsinh_R_dln_R*dln_R_drho;
            double darcsinh_R_dpi = darcsinh_R_dln_I*dln_I_dpi + darcsinh_R_dln_R*dln_R_dpi;
            double darcsinh_R_dsig = darcsinh_R_dln_I*dln_I_dsig + darcsinh_R_dln_R*dln_R_dsig;

            double arcsinh_I = std::atan2(ln_I,ln_R);
            double darcsinh_I_dln_R = -ln_I/(pow(ln_I, 2) + pow(ln_R, 2));
            double darcsinh_I_dln_I = ln_R/(pow(ln_I, 2) + pow(ln_R, 2));
            double darcsinh_I_drho = darcsinh_I_dln_I*dln_I_drho + darcsinh_I_dln_R*dln_R_drho;
            double darcsinh_I_dpi = darcsinh_I_dln_I*dln_I_dpi + darcsinh_I_dln_R*dln_R_dpi;
            double darcsinh_I_dsig = darcsinh_I_dln_I*dln_I_dsig + darcsinh_I_dln_R*dln_R_dsig;

            // fd = 1 + 6 * fb * x * arcsinh(x)
            double fd_R =  1.0 + 6.0 * fb * (arcsinh_R*x_R - arcsinh_I*x_I);
            double dfd_R_dx_I = -6*arcsinh_I*fb;
            double dfd_R_dx_R = 6*arcsinh_R*fb;
            double dfd_R_darcsinh_R = 6*fb*x_R;
            double dfd_R_darcsinh_I = -6*fb*x_I;
            double dfd_R_drho = darcsinh_I_drho*dfd_R_darcsinh_I + darcsinh_R_drho*dfd_R_darcsinh_R + dfd_R_dx_I*dx_I_drho + dfd_R_dx_R*dx_R_drho;
            double dfd_R_dpi = darcsinh_I_dpi*dfd_R_darcsinh_I + darcsinh_R_dpi*dfd_R_darcsinh_R + dfd_R_dx_I*dx_I_dpi + dfd_R_dx_R*dx_R_dpi;
            double dfd_R_dsig = darcsinh_I_dsig*dfd_R_darcsinh_I + darcsinh_R_dsig*dfd_R_darcsinh_R + dfd_R_dx_I*dx_I_dsig + dfd_R_dx_R*dx_R_dsig;

            double fd_I = 6.0 * fb * (arcsinh_R * x_I + arcsinh_I * x_R);
            double dfd_I_dx_I = 6*arcsinh_R*fb;
            double dfd_I_dx_R = 6*arcsinh_I*fb;
            double dfd_I_darcsinh_R = 6*fb*x_I;
            double dfd_I_darcsinh_I = 6*fb*x_R;
            double dfd_I_drho = darcsinh_I_drho*dfd_I_darcsinh_I + darcsinh_R_drho*dfd_I_darcsinh_R + dfd_I_dx_I*dx_I_drho + dfd_I_dx_R*dx_R_drho;
            double dfd_I_dpi = darcsinh_I_dpi*dfd_I_darcsinh_I + darcsinh_R_dpi*dfd_I_darcsinh_R + dfd_I_dx_I*dx_I_dpi + dfd_I_dx_R*dx_R_dpi;
            double dfd_I_dsig = darcsinh_I_dsig*dfd_I_darcsinh_I + darcsinh_R_dsig*dfd_I_darcsinh_R + dfd_I_dx_I*dx_I_dsig + dfd_I_dx_R*dx_R_dsig;

            double R_fd = std::sqrt(fd_R*fd_R + fd_I*fd_I);
            double dR_fd_dfd_R = fd_R/sqrt(pow(fd_I, 2) + pow(fd_R, 2));
            double dR_fd_dfd_I = fd_I/sqrt(pow(fd_I, 2) + pow(fd_R, 2));
            double dR_fd_drho = dR_fd_dfd_I*dfd_I_drho + dR_fd_dfd_R*dfd_R_drho;
            double dR_fd_dpi = dR_fd_dfd_I*dfd_I_dpi + dR_fd_dfd_R*dfd_R_dpi;
            double dR_fd_dsig = dR_fd_dfd_I*dfd_I_dsig + dR_fd_dfd_R*dfd_R_dsig;

            double arg_fd = std::atan2(fd_I, fd_R);
            double darg_fd_dfd_R = -fd_I/(pow(fd_I, 2) + pow(fd_R, 2));
            double darg_fd_dfd_I = fd_R/(pow(fd_I, 2) + pow(fd_R, 2));
            double darg_fd_drho = darg_fd_dfd_I*dfd_I_drho + darg_fd_dfd_R*dfd_R_drho;
            double darg_fd_dpi = darg_fd_dfd_I*dfd_I_dpi + darg_fd_dfd_R*dfd_R_dpi;
            double darg_fd_dsig = darg_fd_dfd_I*dfd_I_dsig + darg_fd_dfd_R*dfd_R_dsig;

            // fzeta = rho**4/3 * x**2 / fd
            fzeta = 2.0 * R * Rx2 / R_fd * std::cos(argr + argx2 - arg_fd);
            double dfzeta_dargx2 = -2*R*Rx2*sin(-arg_fd + argr + argx2)/R_fd;
            double dfzeta_dRx2 = 2*R*cos(-arg_fd + argr + argx2)/R_fd;
            double dfzeta_dargr = -2*R*Rx2*sin(-arg_fd + argr + argx2)/R_fd;
            double dfzeta_dR_fd = -2*R*Rx2*cos(-arg_fd + argr + argx2)/pow(R_fd, 2);
            double dfzeta_dR = 2*Rx2*cos(-arg_fd + argr + argx2)/R_fd;
            double dfzeta_darg_fd = 2*R*Rx2*sin(-arg_fd + argr + argx2)/R_fd;
            dfzeta_drho = dR_drho*dfzeta_dR + dR_fd_drho*dfzeta_dR_fd + dRx2_drho*dfzeta_dRx2 + darg_fd_drho*dfzeta_darg_fd + dargr_drho*dfzeta_dargr + dargx2_drho*dfzeta_dargx2;
            dfzeta_dpi = dR_dpi*dfzeta_dR + dR_fd_dpi*dfzeta_dR_fd + dRx2_dpi*dfzeta_dRx2 + darg_fd_dpi*dfzeta_darg_fd + dargr_dpi*dfzeta_dargr + dargx2_dpi*dfzeta_dargx2;
            dfzeta_dsig = dR_fd_dsig*dfzeta_dR_fd + dRx2_dsig*dfzeta_dRx2 + darg_fd_dsig*dfzeta_darg_fd + dargx2_dsig*dfzeta_dargx2;
        }

        exc[g] = -fb * fzeta/density;

        vrho[2 * g + 0] =  -fb*dfzeta_drho;

        vrho[2 * g + 1] =  -fb*dfzeta_dpi;

        vsigma[3 * g + 0] =  -fb*dfzeta_dsig;
    }
}

}  // namespace pdftb88
