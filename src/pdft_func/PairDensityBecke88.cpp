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

   // Subroutine generated by xc_write in MultiPsi, copyright M.G. Delcey, 2024

   double b88_beta = 0.0041999999999999997;

   double b88_gamma = 6.0;

   for (int32_t g = 0; g < np; g++)
   {
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

      double xt = 0;

      double df_dsig = 0;

      double df_drho = 0;

      double dxt_dsig = 0;

      double dxt_drho = 0;

      double f = 0;

      double df_dpi = 0;

      if (pair_density < -9.9999999999999998e-17)
      {
         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         double zeta = M_SQRT2*sqrt(-pair_density)/density;

         double dzeta_dpi = (1.0/2.0)*M_SQRT2*sqrt(-pair_density)/(density*pair_density);

         double dzeta_drho = -M_SQRT2*sqrt(-pair_density)/pow(density, 2);

         double xa = cbrt(2)*xt/cbrt(zeta + 1);

         double dxa_dxt = cbrt(2)/cbrt(zeta + 1);

         double dxa_dzeta = -1.0/3.0*cbrt(2)*xt/pow(zeta + 1, 4.0/3.0);

         double dxa_dpi = dxa_dzeta*dzeta_dpi;

         double dxa_drho = dxa_dxt*dxt_drho + dxa_dzeta*dzeta_drho;

         double dxa_dsig = dxa_dxt*dxt_dsig;

         double denom_a = b88_beta*b88_gamma*xa*asinh(xa) + 1;

         double ddenom_a_dxa = b88_beta*b88_gamma*xa/sqrt(pow(xa, 2) + 1) + b88_beta*b88_gamma*asinh(xa);

         double ddenom_a_dpi = ddenom_a_dxa*dxa_dpi;

         double ddenom_a_drho = ddenom_a_dxa*dxa_drho;

         double ddenom_a_dsig = ddenom_a_dxa*dxa_dsig;

         double fa = pow(zeta + 1, 2.0/3.0)/denom_a;

         double dfa_ddenom_a = -pow(zeta + 1, 2.0/3.0)/pow(denom_a, 2);

         double dfa_dzeta = (2.0/3.0)/(denom_a*cbrt(zeta + 1));

         double dfa_dsig = ddenom_a_dsig*dfa_ddenom_a;

         double dfa_dpi = ddenom_a_dpi*dfa_ddenom_a + dfa_dzeta*dzeta_dpi;

         double dfa_drho = ddenom_a_drho*dfa_ddenom_a + dfa_dzeta*dzeta_drho;

         double xb = cbrt(2)*xt/cbrt(1 - zeta);

         double dxb_dxt = cbrt(2)/cbrt(1 - zeta);

         double dxb_dzeta = (1.0/3.0)*cbrt(2)*xt/pow(1 - zeta, 4.0/3.0);

         double dxb_dpi = dxb_dzeta*dzeta_dpi;

         double dxb_drho = dxb_dxt*dxt_drho + dxb_dzeta*dzeta_drho;

         double dxb_dsig = dxb_dxt*dxt_dsig;

         double denom_b = b88_beta*b88_gamma*xb*asinh(xb) + 1;

         double ddenom_b_dxb = b88_beta*b88_gamma*xb/sqrt(pow(xb, 2) + 1) + b88_beta*b88_gamma*asinh(xb);

         double ddenom_b_dpi = ddenom_b_dxb*dxb_dpi;

         double ddenom_b_drho = ddenom_b_dxb*dxb_drho;

         double ddenom_b_dsig = ddenom_b_dxb*dxb_dsig;

         double fb = pow(1 - zeta, 2.0/3.0)/denom_b;

         double dfb_ddenom_b = -pow(1 - zeta, 2.0/3.0)/pow(denom_b, 2);

         double dfb_dzeta = -(2.0/3.0)/(denom_b*cbrt(1 - zeta));

         double dfb_dsig = ddenom_b_dsig*dfb_ddenom_b;

         double dfb_dpi = ddenom_b_dpi*dfb_ddenom_b + dfb_dzeta*dzeta_dpi;

         double dfb_drho = ddenom_b_drho*dfb_ddenom_b + dfb_dzeta*dzeta_drho;

         f = fa + fb;

         double df_dfb = 1;

         double df_dfa = 1;

         df_dsig = df_dfa*dfa_dsig + df_dfb*dfb_dsig;

         df_dpi = df_dfa*dfa_dpi + df_dfb*dfb_dpi;

         df_drho = df_dfa*dfa_drho + df_dfb*dfb_drho;

      }
      else if (pair_density < 9.9999999999999998e-17)
      {
         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         double xs = cbrt(2)*xt;

         double dxs_dxt = cbrt(2);

         double dxs_drho = dxs_dxt*dxt_drho;

         double dxs_dsig = dxs_dxt*dxt_dsig;

         double denoms = b88_beta*b88_gamma*xs*asinh(xs) + 1;

         double ddenoms_dxs = b88_beta*b88_gamma*xs/sqrt(pow(xs, 2) + 1) + b88_beta*b88_gamma*asinh(xs);

         double ddenoms_drho = ddenoms_dxs*dxs_drho;

         double ddenoms_dsig = ddenoms_dxs*dxs_dsig;

         double xs2 = pow(xs, 2);

         double dxs2_dxs = 2*xs;

         double dxs2_drho = dxs2_dxs*dxs_drho;

         double dxs2_dsig = dxs2_dxs*dxs_dsig;

         double sqrtx = sqrt(xs2 + 1);

         double dsqrtx_dxs2 = (1.0/2.0)/sqrt(xs2 + 1);

         double dsqrtx_drho = dsqrtx_dxs2*dxs2_drho;

         double dsqrtx_dsig = dsqrtx_dxs2*dxs2_dsig;

         double z2 = (1.0/9.0)*(-2*b88_beta*b88_gamma*xs2/sqrtx + b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 3) - 2*denoms + 2*pow(b88_beta*b88_gamma*xs2/sqrtx + denoms - 1, 2)/denoms)/pow(denoms, 2);

         double dz2_ddenoms = (1.0/9.0)*(-2 + 2*(2*b88_beta*b88_gamma*xs2/sqrtx + 2*denoms - 2)/denoms - 2*pow(b88_beta*b88_gamma*xs2/sqrtx + denoms - 1, 2)/pow(denoms, 2))/pow(denoms, 2) - 2.0/9.0*(-2*b88_beta*b88_gamma*xs2/sqrtx + b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 3) - 2*denoms + 2*pow(b88_beta*b88_gamma*xs2/sqrtx + denoms - 1, 2)/denoms)/pow(denoms, 3);

         double dz2_dxs = (4.0/9.0)*b88_beta*b88_gamma*pow(xs, 3)/(pow(denoms, 2)*pow(sqrtx, 3));

         double dz2_dsqrtx = (1.0/9.0)*(2*b88_beta*b88_gamma*xs2/pow(sqrtx, 2) - 3*b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 4) - 4*b88_beta*b88_gamma*xs2*(b88_beta*b88_gamma*xs2/sqrtx + denoms - 1)/(denoms*pow(sqrtx, 2)))/pow(denoms, 2);

         double dz2_dxs2 = (1.0/9.0)*(-2*b88_beta*b88_gamma/sqrtx + 4*b88_beta*b88_gamma*(b88_beta*b88_gamma*xs2/sqrtx + denoms - 1)/(denoms*sqrtx))/pow(denoms, 2);

         double dz2_drho = ddenoms_drho*dz2_ddenoms + dsqrtx_drho*dz2_dsqrtx + dxs2_drho*dz2_dxs2 + dxs_drho*dz2_dxs;

         double dz2_dsig = ddenoms_dsig*dz2_ddenoms + dsqrtx_dsig*dz2_dsqrtx + dxs2_dsig*dz2_dxs2 + dxs_dsig*dz2_dxs;

         double zeta2 = -2*pair_density/pow(density, 2);

         double dzeta2_dpi = -2/pow(density, 2);

         double dzeta2_drho = 4*pair_density/pow(density, 3);

         f = z2*zeta2 + 2/denoms;

         double df_dzeta2 = z2;

         double df_ddenoms = -2/pow(denoms, 2);

         double df_dz2 = zeta2;

         df_dpi = df_dzeta2*dzeta2_dpi;

         df_drho = ddenoms_drho*df_ddenoms + df_dz2*dz2_drho + df_dzeta2*dzeta2_drho;

         df_dsig = ddenoms_dsig*df_ddenoms + df_dz2*dz2_dsig;

      }
      else
      {
         double eta = M_SQRT2*sqrt(pair_density)/density;

         double deta_dpi = (1.0/2.0)*M_SQRT2/(density*sqrt(pair_density));

         double deta_drho = -M_SQRT2*sqrt(pair_density)/pow(density, 2);

         double theta = -1.0/3.0*atan(eta);

         double dtheta_deta = -(1.0/3.0)/(pow(eta, 2) + 1);

         double dtheta_dpi = deta_dpi*dtheta_deta;

         double dtheta_drho = deta_drho*dtheta_deta;

         double r = pow(pow(eta, 2) + 1, -1.0/6.0);

         double dr_deta = -1.0/3.0*eta/pow(pow(eta, 2) + 1, 7.0/6.0);

         double dr_dpi = deta_dpi*dr_deta;

         double dr_drho = deta_drho*dr_deta;

         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         double xR = cbrt(2)*r*xt*cos(theta);

         double dxR_dxt = cbrt(2)*r*cos(theta);

         double dxR_dr = cbrt(2)*xt*cos(theta);

         double dxR_dtheta = -cbrt(2)*r*xt*sin(theta);

         double dxR_dsig = dxR_dxt*dxt_dsig;

         double dxR_dpi = dr_dpi*dxR_dr + dtheta_dpi*dxR_dtheta;

         double dxR_drho = dr_drho*dxR_dr + dtheta_drho*dxR_dtheta + dxR_dxt*dxt_drho;

         double xI = cbrt(2)*r*xt*sin(theta);

         double dxI_dxt = cbrt(2)*r*sin(theta);

         double dxI_dr = cbrt(2)*xt*sin(theta);

         double dxI_dtheta = cbrt(2)*r*xt*cos(theta);

         double dxI_dsig = dxI_dxt*dxt_dsig;

         double dxI_dpi = dr_dpi*dxI_dr + dtheta_dpi*dxI_dtheta;

         double dxI_drho = dr_drho*dxI_dr + dtheta_drho*dxI_dtheta + dxI_dxt*dxt_drho;

         //arcsinh(a) = ln(a + sqrt(a^2+1))
         double x2R = -pow(xI, 2) + pow(xR, 2) + 1;

         double dx2R_dxI = -2*xI;

         double dx2R_dxR = 2*xR;

         double dx2R_dsig = dx2R_dxI*dxI_dsig + dx2R_dxR*dxR_dsig;

         double dx2R_dpi = dx2R_dxI*dxI_dpi + dx2R_dxR*dxR_dpi;

         double dx2R_drho = dx2R_dxI*dxI_drho + dx2R_dxR*dxR_drho;

         double x2I = 2*xI*xR;

         double dx2I_dxI = 2*xR;

         double dx2I_dxR = 2*xI;

         double dx2I_dsig = dx2I_dxI*dxI_dsig + dx2I_dxR*dxR_dsig;

         double dx2I_dpi = dx2I_dxI*dxI_dpi + dx2I_dxR*dxR_dpi;

         double dx2I_drho = dx2I_dxI*dxI_drho + dx2I_dxR*dxR_drho;

         double rsqrt = pow(pow(x2I, 2) + pow(x2R, 2), 1.0/4.0);

         double drsqrt_dx2I = (1.0/2.0)*x2I/pow(pow(x2I, 2) + pow(x2R, 2), 3.0/4.0);

         double drsqrt_dx2R = (1.0/2.0)*x2R/pow(pow(x2I, 2) + pow(x2R, 2), 3.0/4.0);

         double drsqrt_dsig = drsqrt_dx2I*dx2I_dsig + drsqrt_dx2R*dx2R_dsig;

         double drsqrt_dpi = drsqrt_dx2I*dx2I_dpi + drsqrt_dx2R*dx2R_dpi;

         double drsqrt_drho = drsqrt_dx2I*dx2I_drho + drsqrt_dx2R*dx2R_drho;

         double argsqrt = 0.5*atan2(x2I, x2R);

         double dargsqrt_dx2I = 0.5*x2R/(pow(x2I, 2) + pow(x2R, 2));

         double dargsqrt_dx2R = -0.5*x2I/(pow(x2I, 2) + pow(x2R, 2));

         double dargsqrt_dsig = dargsqrt_dx2I*dx2I_dsig + dargsqrt_dx2R*dx2R_dsig;

         double dargsqrt_dpi = dargsqrt_dx2I*dx2I_dpi + dargsqrt_dx2R*dx2R_dpi;

         double dargsqrt_drho = dargsqrt_dx2I*dx2I_drho + dargsqrt_dx2R*dx2R_drho;

         double lnR = rsqrt*cos(argsqrt) + xR;

         double dlnR_drsqrt = cos(argsqrt);

         double dlnR_dargsqrt = -rsqrt*sin(argsqrt);

         double dlnR_dxR = 1;

         double dlnR_dsig = dargsqrt_dsig*dlnR_dargsqrt + dlnR_drsqrt*drsqrt_dsig + dlnR_dxR*dxR_dsig;

         double dlnR_dpi = dargsqrt_dpi*dlnR_dargsqrt + dlnR_drsqrt*drsqrt_dpi + dlnR_dxR*dxR_dpi;

         double dlnR_drho = dargsqrt_drho*dlnR_dargsqrt + dlnR_drsqrt*drsqrt_drho + dlnR_dxR*dxR_drho;

         double lnI = rsqrt*sin(argsqrt) + xI;

         double dlnI_dargsqrt = rsqrt*cos(argsqrt);

         double dlnI_dxI = 1;

         double dlnI_drsqrt = sin(argsqrt);

         double dlnI_dsig = dargsqrt_dsig*dlnI_dargsqrt + dlnI_drsqrt*drsqrt_dsig + dlnI_dxI*dxI_dsig;

         double dlnI_dpi = dargsqrt_dpi*dlnI_dargsqrt + dlnI_drsqrt*drsqrt_dpi + dlnI_dxI*dxI_dpi;

         double dlnI_drho = dargsqrt_drho*dlnI_dargsqrt + dlnI_drsqrt*drsqrt_drho + dlnI_dxI*dxI_drho;

         double arcsinR = 0.5*log(pow(lnI, 2) + pow(lnR, 2));

         double darcsinR_dlnI = 1.0*lnI/(pow(lnI, 2) + pow(lnR, 2));

         double darcsinR_dlnR = 1.0*lnR/(pow(lnI, 2) + pow(lnR, 2));

         double darcsinR_dsig = darcsinR_dlnI*dlnI_dsig + darcsinR_dlnR*dlnR_dsig;

         double darcsinR_dpi = darcsinR_dlnI*dlnI_dpi + darcsinR_dlnR*dlnR_dpi;

         double darcsinR_drho = darcsinR_dlnI*dlnI_drho + darcsinR_dlnR*dlnR_drho;

         double arcsinI = atan2(lnI, lnR);

         double darcsinI_dlnI = lnR/(pow(lnI, 2) + pow(lnR, 2));

         double darcsinI_dlnR = -lnI/(pow(lnI, 2) + pow(lnR, 2));

         double darcsinI_dsig = darcsinI_dlnI*dlnI_dsig + darcsinI_dlnR*dlnR_dsig;

         double darcsinI_dpi = darcsinI_dlnI*dlnI_dpi + darcsinI_dlnR*dlnR_dpi;

         double darcsinI_drho = darcsinI_dlnI*dlnI_drho + darcsinI_dlnR*dlnR_drho;

         double demomR = b88_beta*b88_gamma*(-arcsinI*xI + arcsinR*xR) + 1;

         double ddemomR_dxR = arcsinR*b88_beta*b88_gamma;

         double ddemomR_darcsinR = b88_beta*b88_gamma*xR;

         double ddemomR_dxI = -arcsinI*b88_beta*b88_gamma;

         double ddemomR_darcsinI = -b88_beta*b88_gamma*xI;

         double ddemomR_dsig = darcsinI_dsig*ddemomR_darcsinI + darcsinR_dsig*ddemomR_darcsinR + ddemomR_dxI*dxI_dsig + ddemomR_dxR*dxR_dsig;

         double ddemomR_dpi = darcsinI_dpi*ddemomR_darcsinI + darcsinR_dpi*ddemomR_darcsinR + ddemomR_dxI*dxI_dpi + ddemomR_dxR*dxR_dpi;

         double ddemomR_drho = darcsinI_drho*ddemomR_darcsinI + darcsinR_drho*ddemomR_darcsinR + ddemomR_dxI*dxI_drho + ddemomR_dxR*dxR_drho;

         double demomI = b88_beta*b88_gamma*(arcsinI*xR + arcsinR*xI);

         double ddemomI_darcsinR = b88_beta*b88_gamma*xI;

         double ddemomI_dxR = arcsinI*b88_beta*b88_gamma;

         double ddemomI_darcsinI = b88_beta*b88_gamma*xR;

         double ddemomI_dxI = arcsinR*b88_beta*b88_gamma;

         double ddemomI_dsig = darcsinI_dsig*ddemomI_darcsinI + darcsinR_dsig*ddemomI_darcsinR + ddemomI_dxI*dxI_dsig + ddemomI_dxR*dxR_dsig;

         double ddemomI_dpi = darcsinI_dpi*ddemomI_darcsinI + darcsinR_dpi*ddemomI_darcsinR + ddemomI_dxI*dxI_dpi + ddemomI_dxR*dxR_dpi;

         double ddemomI_drho = darcsinI_drho*ddemomI_darcsinI + darcsinR_drho*ddemomI_darcsinR + ddemomI_dxI*dxI_drho + ddemomI_dxR*dxR_drho;

         double rdenom = sqrt(pow(demomI, 2) + pow(demomR, 2));

         double drdenom_ddemomR = demomR/sqrt(pow(demomI, 2) + pow(demomR, 2));

         double drdenom_ddemomI = demomI/sqrt(pow(demomI, 2) + pow(demomR, 2));

         double drdenom_dsig = ddemomI_dsig*drdenom_ddemomI + ddemomR_dsig*drdenom_ddemomR;

         double drdenom_dpi = ddemomI_dpi*drdenom_ddemomI + ddemomR_dpi*drdenom_ddemomR;

         double drdenom_drho = ddemomI_drho*drdenom_ddemomI + ddemomR_drho*drdenom_ddemomR;

         double argdenom = atan2(demomI, demomR);

         double dargdenom_ddemomR = -demomI/(pow(demomI, 2) + pow(demomR, 2));

         double dargdenom_ddemomI = demomR/(pow(demomI, 2) + pow(demomR, 2));

         double dargdenom_dsig = dargdenom_ddemomI*ddemomI_dsig + dargdenom_ddemomR*ddemomR_dsig;

         double dargdenom_dpi = dargdenom_ddemomI*ddemomI_dpi + dargdenom_ddemomR*ddemomR_dpi;

         double dargdenom_drho = dargdenom_ddemomI*ddemomI_drho + dargdenom_ddemomR*ddemomR_drho;

         f = 2*cos(argdenom + 2*theta)/(pow(r, 2)*rdenom);

         double df_dr = -4*cos(argdenom + 2*theta)/(pow(r, 3)*rdenom);

         double df_dtheta = -4*sin(argdenom + 2*theta)/(pow(r, 2)*rdenom);

         double df_dargdenom = -2*sin(argdenom + 2*theta)/(pow(r, 2)*rdenom);

         double df_drdenom = -2*cos(argdenom + 2*theta)/(pow(r, 2)*pow(rdenom, 2));

         df_dsig = dargdenom_dsig*df_dargdenom + df_drdenom*drdenom_dsig;

         df_dpi = dargdenom_dpi*df_dargdenom + df_dr*dr_dpi + df_drdenom*drdenom_dpi + df_dtheta*dtheta_dpi;

         df_drho = dargdenom_drho*df_dargdenom + df_dr*dr_drho + df_drdenom*drdenom_drho + df_dtheta*dtheta_drho;

      }
      double xt2 = pow(density, -2.6666666666666665)*sig;

      double dxt2_drho = -2.6666666666666665*pow(density, -3.6666666666666665)*sig;

      double dxt2_dsig = pow(density, -2.6666666666666665);

      exc[g] =  -1.0/2.0*cbrt(2)*b88_beta*cbrt(density)*f*xt2;

      double dExc_dxt2 = -1.0/2.0*cbrt(2)*b88_beta*pow(density, 4.0/3.0)*f;

      double dExc_df = -1.0/2.0*cbrt(2)*b88_beta*pow(density, 4.0/3.0)*xt2;

      vrho[2 * g + 1] =  dExc_df*df_dpi;

      vrho[2 * g + 0] =  -2.0/3.0*cbrt(2)*b88_beta*cbrt(density)*f*xt2 + dExc_df*df_drho + dExc_dxt2*dxt2_drho;

      vsigma[3 * g + 0] =  dExc_df*df_dsig + dExc_dxt2*dxt2_dsig;

      //Currently, no explicit dependence
      vsigma[3 * g + 1] = 0.0;
      vsigma[3 * g + 2] = 0.0;

   }
}

}  // namespace pdftb88
