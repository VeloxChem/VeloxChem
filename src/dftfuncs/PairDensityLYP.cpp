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

#include "PairDensityLYP.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftlyp {  // pdftlyp namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
{
   //Constants
   double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
   double cf = (3.0/10.0)*pow(3, 2.0/3.0)*pow(M_PI, 4.0/3.0);

   double f83 = 8.0/3.0;

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
         vsigma[2 * g + 0] = 0.0;
         continue;
      }

      double pair_density = rho[2 * g + 1];

      double sig = sigma[3 * g + 0];

      double rho2 = pow(density, 2);
      double drho2_drho = 2*density;

      double rhom13 = pow(density, -1.0/3.0);
      double drhom13_drho = -(1.0/3.0)/pow(density, 4.0/3.0);

      double denom2 = d*rhom13 + 1;
      double ddenom2_drho = d*drhom13_drho;

      double omega = exp(-c*rhom13)/(denom2*pow(density, 11.0/3.0));
      double domega_drho = -c*drhom13_drho*exp(-c*rhom13)/(denom2*pow(density, 11.0/3.0)) - ddenom2_drho*exp(-c*rhom13)/(pow(denom2, 2)*pow(density, 11.0/3.0)) - 11.0/3.0*exp(-c*rhom13)/(denom2*pow(density, 14.0/3.0));

      double delta = rhom13*(c + d/denom2);
      double ddelta_drho = -d*ddenom2_drho*rhom13/pow(denom2, 2) + drhom13_drho*(c + d/denom2);

      double t2 = (1.0/18.0)*sig*(47 - 7*delta);
      double dt2_drho = -7.0/18.0*ddelta_drho*sig;
      double dt2_dsig = 47.0/18.0 - 7.0/18.0*delta;

      double t5 = -2.0/3.0*rho2*sig;
      double dt5_drho = -2.0/3.0*drho2_drho*sig;
      double dt5_dsig = -2.0/3.0*rho2;

      double delta2 = -2*pair_density;
      double ddelta2_dpi = -2;

      double zeta2 = delta2/rho2;
      double dzeta2_drho = -delta2*drho2_drho/pow(rho2, 2);
      double dzeta2_dpi = ddelta2_dpi/rho2;

      double gradR = (1.0/4.0)*sig*(zeta2 + 1);
      double dgradR_dzeta2 = (1.0/4.0)*sig;
      double dgradR_drho = dgradR_dzeta2*dzeta2_drho;
      double dgradR_dpi = dgradR_dzeta2*dzeta2_dpi;
      double dgradR_dsig = (1.0/4.0)*zeta2 + 1.0/4.0;

      double t3 = 2*gradR*((1.0/18.0)*delta - 5.0/2.0);
      double dt3_dgradR = (1.0/9.0)*delta - 5;
      double dt3_drho = (1.0/9.0)*ddelta_drho*gradR + dgradR_drho*dt3_dgradR;
      double dt3_dpi = dgradR_dpi*dt3_dgradR;
      double dt3_dsig = dgradR_dsig*dt3_dgradR;

      double sA = gradR + 0.5 * sig*zeta2;
      double dsA_drho = dgradR_drho + 0.5 * sig * dzeta2_drho;
      double dsA_dpi = dgradR_dpi + 0.5 * sig * dzeta2_dpi;
      double dsA_dsig = dgradR_dsig + 0.5 * zeta2;

      double t1 = 0;
      double dt1_drho = 0;
      double dt1_dpi = 0;

      double t6 = 0;
      double dt6_drho = 0;
      double dt6_dpi = 0;
      double dt6_dsig = 0;

      // Real part
      if (pair_density < -1.0e-30)
      {
         double zeta = sqrt(-2.0*pair_density)/density;
         double dzeta_drho = -sqrt(-2.0*pair_density)/pow(density, 2);
         double dzeta_dpi = (1.0/2.0)*sqrt(-2.0*pair_density)/(density*pair_density);

         t1 = 2*cf*pow(density, f83)*(pow(1 - zeta, f83) + pow(zeta + 1, f83));
         double dt1_dzeta = 2*cf*pow(density, f83)*(-f83*pow(1 - zeta, 5.0/3.0) + f83*pow(zeta + 1, 5.0/3.0));
         dt1_drho = (16.0/3.0)*cf*pow(density, 5.0/3.0)*(pow(1 - zeta, f83) + pow(zeta + 1, f83)) + dt1_dzeta*dzeta_drho;
         dt1_dpi = dt1_dzeta*dzeta_dpi;

         t6 = (gradR - 1.0/2.0*sig*zeta)*(-1.0/4.0*rho2*(2*zeta + zeta2 + 1) + (2.0/3.0)*rho2) + (gradR + (1.0/2.0)*sig*zeta)*(-1.0/4.0*rho2*(-2*zeta + zeta2 + 1) + (2.0/3.0)*rho2);
         double dt6_dgradR = -1.0/4.0*rho2*(-2*zeta + zeta2 + 1) - 1.0/4.0*rho2*(2*zeta + zeta2 + 1) + (4.0/3.0)*rho2;
         double dt6_dzeta = -1.0/2.0*rho2*(gradR - 1.0/2.0*sig*zeta) + (1.0/2.0)*rho2*(gradR + (1.0/2.0)*sig*zeta) + (1.0/2.0)*sig*(-1.0/4.0*rho2*(-2*zeta + zeta2 + 1) + (2.0/3.0)*rho2) - 1.0/2.0*sig*(-1.0/4.0*rho2*(2*zeta + zeta2 + 1) + (2.0/3.0)*rho2);
         double dt6_dzeta2 = -1.0/4.0*rho2*(gradR - 1.0/2.0*sig*zeta) - 1.0/4.0*rho2*(gradR + (1.0/2.0)*sig*zeta);
         dt6_drho = dgradR_drho*dt6_dgradR + drho2_drho*((gradR - 1.0/2.0*sig*zeta)*(-1.0/2.0*zeta - 1.0/4.0*zeta2 + 5.0/12.0) + (gradR + (1.0/2.0)*sig*zeta)*((1.0/2.0)*zeta - 1.0/4.0*zeta2 + 5.0/12.0)) + dt6_dzeta*dzeta_drho + dt6_dzeta2*dzeta2_drho;
         dt6_dpi = dgradR_dpi*dt6_dgradR + dt6_dzeta*dzeta_dpi + dt6_dzeta2*dzeta2_dpi;
         dt6_dsig = dgradR_dsig*dt6_dgradR + (1.0/2.0)*zeta*(-1.0/4.0*rho2*(-2*zeta + zeta2 + 1) + (2.0/3.0)*rho2) - 1.0/2.0*zeta*(-1.0/4.0*rho2*(2*zeta + zeta2 + 1) + (2.0/3.0)*rho2);
      }
      // pi=0 special case
      else if (pair_density < 1.0e-30)
      {
         t1 = 4*cf*pow(density, f83);
         dt1_dpi = -160.0/9.0*cf*pow(density, f83)/pow(density, 2);
         dt1_drho = (32.0/3.0)*cf*pow(density, 5.0/3.0);

         t6 = 2*gradR*rho2 * (-1.0/4.0*(zeta2 + 1) + (2.0/3.0));
         dt6_dpi = -sig  + 0.25*sig * (zeta2 - 5.0/3.0) + gradR;
         double dt6_dgradR = -1.0/2.0*rho2*(zeta2 + 1) + (4.0/3.0)*rho2;
         double dt6_dzeta2 = -1.0/2.0*gradR*rho2;
         dt6_drho = dgradR_drho*dt6_dgradR + 2*drho2_drho*gradR*(5.0/12.0 - 1.0/4.0*zeta2) + dt6_dzeta2*dzeta2_drho;
         dt6_dsig = dgradR_dsig*dt6_dgradR;
      }
      // imaginary part
      else
      {
         double eta = sqrt(2.0*pair_density)/density;
         double deta_drho = -sqrt(2.0*pair_density)/pow(density, 2);
         double deta_dpi = 1.0/(density*sqrt(2.0*pair_density));

         t1 = 4*cf*pow(1.0*pow(eta, 2) + 1, 4.0/3.0)*cos(f83*atan(eta))*pow(density, f83);
         double eta_arg = f83*atan(eta);
         double eta_norm = cbrt(1 + pow(eta, 2));
         double dt1_deta = 32.0/3.0*cf*eta*eta_norm*cos(eta_arg)*pow(density, f83) - 32.0/3.0*cf*cbrt(pow(eta, 2) + 1)*sin(eta_arg)*pow(density, f83);
         dt1_drho = (32.0/3.0)*cf*pow(pow(eta, 2) + 1, 4.0/3.0)*cos(eta_arg)*pow(density, 5.0/3.0) + deta_drho*dt1_deta;
         dt1_dpi = deta_dpi*dt1_deta;

         t6 = -0.5*pow(eta, 2)*rho2*sig - 1.0/2.0*gradR*rho2*(zeta2 + 1) + (4.0/3.0)*gradR*rho2;
         double dt6_deta = - eta*rho2*sig;
         double dt6_dgradR = -1.0/2.0*rho2*(zeta2 + 1) + (4.0/3.0)*rho2;
         double dt6_dzeta2 = -1.0/2.0*gradR*rho2;
         dt6_drho = deta_drho*dt6_deta + dgradR_drho*dt6_dgradR + drho2_drho*(-0.5*pow(eta, 2)*sig - 1.0/2.0*gradR*(zeta2 + 1) + (4.0/3.0)*gradR) + dt6_dzeta2*dzeta2_drho;
         dt6_dpi = deta_dpi*dt6_deta + dgradR_dpi*dt6_dgradR + dt6_dzeta2*dzeta2_dpi;
         dt6_dsig = dgradR_dsig*dt6_dgradR - 0.5*pow(eta, 2)*rho2;
      }

      double t4 = sA*(11.0/9.0 - 1.0/9.0*delta);
      double dt4_dsA = 11.0/9.0 - 1.0/9.0*delta;
      double dt4_drho = -1.0/9.0*ddelta_drho*sA + dsA_drho*dt4_dsA;
      double dt4_dpi = dsA_dpi*dt4_dsA;
      double dt4_dsig = dsA_dsig*dt4_dsA;

      double rarb = (1.0/4.0)*pow(density, 2)*(-delta2/rho2 + 1);
      double drarb_drho = (1.0/4.0)*delta2*pow(density, 2)*drho2_drho/pow(rho2, 2) + (1.0/2.0)*density*(-delta2/rho2 + 1);
      double drarb_dpi = -1.0/4.0*ddelta2_dpi*pow(density, 2)/rho2;

      exc[g] =  -a*(b*omega*(rarb*(t1 + t2 + t3 + t4) + t5 + t6) + 4*rarb/(denom2*density))/density;
      double dExc_dt4 = -a*b*omega*rarb;
      double dExc_dt3 = -a*b*omega*rarb;
      double dExc_dt5 = -a*b*omega;
      double dExc_dt2 = -a*b*omega*rarb;
      double dExc_drarb = -a*(b*omega*(t1 + t2 + t3 + t4) + 4/(denom2*density));
      double dExc_dt6 = -a*b*omega;
      double dExc_dt1 = -a*b*omega*rarb;

      vrho[2 * g + 0] =  -a*b*domega_drho*(rarb*(t1 + t2 + t3 + t4) + t5 + t6) + 4*a*ddenom2_drho*rarb/(pow(denom2, 2)*density) + 4*a*rarb/(denom2*pow(density, 2)) + dExc_drarb*drarb_drho + dExc_dt1*dt1_drho + dExc_dt2*dt2_drho + dExc_dt3*dt3_drho + dExc_dt4*dt4_drho + dExc_dt5*dt5_drho + dExc_dt6*dt6_drho;
      vrho[2 * g + 1] =  dExc_drarb*drarb_dpi + dExc_dt1*dt1_dpi + dExc_dt3*dt3_dpi + dExc_dt4*dt4_dpi + dExc_dt6*dt6_dpi;

      vsigma[3 * g + 0] =  dExc_dt2*dt2_dsig + dExc_dt3*dt3_dsig + dExc_dt4*dt4_dsig + dExc_dt5*dt5_dsig + dExc_dt6*dt6_dsig;
   }
}

}  // namespace pdftpbe_c
