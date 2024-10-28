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

#include "PairDensityB88_erf.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftb88_erf {  // pdftb88_erf namespace

void
compute_exc_vxc(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma, const double mu)
{

   // Subroutine generated by xc_write in MultiPsi, copyright M.G. Delcey, 2024

   double x_factor_c = (3.0/4.0)*cbrt(6)/cbrt(M_PI);

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

      double datt_erf0_drho = 0;

      double db88_denom_drho = 0;

      double f = 0;

      double xt = 0;

      double xs = 0;

      double db88_denom_dsig = 0;

      double dxs_dsig = 0;

      double dityh_f_drho = 0;

      double datt_erf_aux3_drho = 0;

      double dxt_drho = 0;

      double datt_erf0_datt_erf_aux3 = 0;

      double b88_f = 0;

      double df_dpi = 0;

      double df_drho = 0;

      double dityh_f_dsig = 0;

      double dityh_k_GGA_dsig = 0;

      double datt_erf_aux1_dityh_k_GGA = 0;

      double dityh_k_GGA_drho = 0;

      double db88_denom_dxs = 0;

      double att_erf_aux1 = 0;

      double dityh_f_db88_denom = 0;

      double datt_erf0_datt_erf_aux1 = 0;

      double db88_f_db88_denom = 0;

      double dxs_drho = 0;

      double ityh_f = 0;

      double db88_f_dxs = 0;

      double datt_erf_aux3_dityh_k_GGA = 0;

      double drho_s_drho = 0;

      double datt_erf_aux1_drho = 0;

      double datt_erf_aux1_dsig = 0;

      double att_erf_aux3 = 0;

      double dxt_dsig = 0;

      double datt_erf_aux3_dsig = 0;

      double datt_erf_aux3_datt_erf_aux2 = 0;

      double datt_erf0_dsig = 0;

      double datt_erf0_dityh_k_GGA = 0;

      double datt_erf0_datt_erf_aux2 = 0;

      double df_dsig = 0;

      double b88_denom = 0;

      double att_erf0 = 0;

      double rho_s = 0;

      double dxs_dxt = 0;

      double db88_f_drho = 0;

      double dityh_k_GGA_db88_f = 0;

      double dityh_f_datt_erf0 = 0;

      double db88_f_dsig = 0;

      double ityh_k_GGA = 0;

      if (fabs(pair_density) > 1.0000000000000001e-32)
      {
         double dzeta_abs_dpi = 0;

         double zeta_abs = 0;

         double dzeta_abs_drho = 0;

         if (pair_density < 0)
         {
            double zeta = M_SQRT2*sqrt(-pair_density)/density;

            double dzeta_drho = -M_SQRT2*sqrt(-pair_density)/pow(density, 2);

            double dzeta_dpi = (1.0/2.0)*M_SQRT2*sqrt(-pair_density)/(density*pair_density);

            zeta_abs = zeta;

            double dzeta_abs_dzeta = 1;

            dzeta_abs_drho = dzeta_abs_dzeta*dzeta_drho;

            dzeta_abs_dpi = dzeta_abs_dzeta*dzeta_dpi;

         }
         else
         {
            double eta = M_SQRT2*sqrt(pair_density)/density;

            double deta_drho = -M_SQRT2*sqrt(pair_density)/pow(density, 2);

            double deta_dpi = (1.0/2.0)*M_SQRT2/(density*sqrt(pair_density));

            zeta_abs = eta;

            double dzeta_abs_deta = 1;

            dzeta_abs_drho = deta_drho*dzeta_abs_deta;

            dzeta_abs_dpi = deta_dpi*dzeta_abs_deta;

         }
         rho_s = 0.5*density*(zeta_abs + 1);

         double drho_s_dzeta_abs = 0.5*density;

         drho_s_drho = drho_s_dzeta_abs*dzeta_abs_drho + 0.5*zeta_abs + 0.5;

         double drho_s_dpi = drho_s_dzeta_abs*dzeta_abs_dpi;

         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         xs = cbrt(2)*xt/cbrt(zeta_abs + 1);

         dxs_dxt = cbrt(2)/cbrt(zeta_abs + 1);

         double dxs_dzeta_abs = -1.0/3.0*cbrt(2)*xt/pow(zeta_abs + 1, 4.0/3.0);

         dxs_drho = dxs_dxt*dxt_drho + dxs_dzeta_abs*dzeta_abs_drho;

         dxs_dsig = dxs_dxt*dxt_dsig;

         double dxs_dpi = dxs_dzeta_abs*dzeta_abs_dpi;

         b88_denom = b88_beta*b88_gamma*xs*asinh(xs) + 1;

         db88_denom_dxs = b88_beta*b88_gamma*xs/sqrt(pow(xs, 2) + 1) + b88_beta*b88_gamma*asinh(xs);

         db88_denom_drho = db88_denom_dxs*dxs_drho;

         db88_denom_dsig = db88_denom_dxs*dxs_dsig;

         double db88_denom_dpi = db88_denom_dxs*dxs_dpi;

         b88_f = b88_beta*pow(xs, 2)/(b88_denom*x_factor_c) + 1;

         db88_f_db88_denom = -b88_beta*pow(xs, 2)/(pow(b88_denom, 2)*x_factor_c);

         db88_f_dxs = 2*b88_beta*xs/(b88_denom*x_factor_c);

         db88_f_drho = db88_denom_drho*db88_f_db88_denom + db88_f_dxs*dxs_drho;

         db88_f_dsig = db88_denom_dsig*db88_f_db88_denom + db88_f_dxs*dxs_dsig;

         double db88_f_dpi = db88_denom_dpi*db88_f_db88_denom + db88_f_dxs*dxs_dpi;

         ityh_k_GGA = (3.0/2.0)*M_SQRT2*sqrt(M_PI)*cbrt(rho_s)*sqrt(1/(b88_f*x_factor_c));

         double dityh_k_GGA_drho_s = (1.0/2.0)*M_SQRT2*sqrt(M_PI)*sqrt(1/(b88_f*x_factor_c))/pow(rho_s, 2.0/3.0);

         dityh_k_GGA_db88_f = -3.0/4.0*M_SQRT2*sqrt(M_PI)*cbrt(rho_s)*sqrt(1/(b88_f*x_factor_c))/b88_f;

         dityh_k_GGA_drho = db88_f_drho*dityh_k_GGA_db88_f + dityh_k_GGA_drho_s*drho_s_drho;

         double dityh_k_GGA_dpi = db88_f_dpi*dityh_k_GGA_db88_f + dityh_k_GGA_drho_s*drho_s_dpi;

         dityh_k_GGA_dsig = db88_f_dsig*dityh_k_GGA_db88_f;

         att_erf_aux1 = sqrt(M_PI)*erf(ityh_k_GGA/mu);

         datt_erf_aux1_dityh_k_GGA = 2*exp(-pow(ityh_k_GGA, 2)/pow(mu, 2))/mu;

         datt_erf_aux1_drho = datt_erf_aux1_dityh_k_GGA*dityh_k_GGA_drho;

         datt_erf_aux1_dsig = datt_erf_aux1_dityh_k_GGA*dityh_k_GGA_dsig;

         double datt_erf_aux1_dpi = datt_erf_aux1_dityh_k_GGA*dityh_k_GGA_dpi;

         double datt_erf_aux2_dsig = 0;

         double datt_erf_aux2_drho = 0;

         double datt_erf_aux2_dityh_k_GGA = 0;

         double att_erf_aux2 = 0;

         double datt_erf_aux2_dpi = 0;

         if ((1.0/2.0)*mu/ityh_k_GGA < 5)
         {
            att_erf_aux2 = -1 + exp(-pow(ityh_k_GGA, 2)/pow(mu, 2));

            datt_erf_aux2_dityh_k_GGA = -2*ityh_k_GGA*exp(-pow(ityh_k_GGA, 2)/pow(mu, 2))/pow(mu, 2);

            datt_erf_aux2_drho = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_drho;

            datt_erf_aux2_dsig = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dsig;

            datt_erf_aux2_dpi = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dpi;

         }
         else
         {
            att_erf_aux2 = -pow(ityh_k_GGA, 2)*(-2*pow(ityh_k_GGA, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_dityh_k_GGA = 4*pow(ityh_k_GGA, 3)/pow(mu, 4) - 2*ityh_k_GGA*(-2*pow(ityh_k_GGA, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_drho = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_drho;

            datt_erf_aux2_dsig = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dsig;

            datt_erf_aux2_dpi = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dpi;

         }
         att_erf_aux3 = (1.0/2.0)*att_erf_aux2*pow(mu, 2)/pow(ityh_k_GGA, 2) + 1.0/2.0;

         datt_erf_aux3_dityh_k_GGA = -att_erf_aux2*pow(mu, 2)/pow(ityh_k_GGA, 3);

         datt_erf_aux3_datt_erf_aux2 = (1.0/2.0)*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf_aux3_drho = datt_erf_aux2_drho*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_dityh_k_GGA*dityh_k_GGA_drho;

         datt_erf_aux3_dsig = datt_erf_aux2_dsig*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_dityh_k_GGA*dityh_k_GGA_dsig;

         double datt_erf_aux3_dpi = datt_erf_aux2_dpi*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_dityh_k_GGA*dityh_k_GGA_dpi;

         att_erf0 = 1 - 4.0/3.0*mu*(att_erf_aux1 + mu*(att_erf_aux2 - att_erf_aux3)/ityh_k_GGA)/ityh_k_GGA;

         datt_erf0_dityh_k_GGA = (4.0/3.0)*mu*(att_erf_aux1 + mu*(att_erf_aux2 - att_erf_aux3)/ityh_k_GGA)/pow(ityh_k_GGA, 2) + (4.0/3.0)*pow(mu, 2)*(att_erf_aux2 - att_erf_aux3)/pow(ityh_k_GGA, 3);

         datt_erf0_datt_erf_aux2 = -4.0/3.0*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf0_datt_erf_aux3 = (4.0/3.0)*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf0_datt_erf_aux1 = -4.0/3.0*mu/ityh_k_GGA;

         datt_erf0_drho = datt_erf0_datt_erf_aux1*datt_erf_aux1_drho + datt_erf0_datt_erf_aux2*datt_erf_aux2_drho + datt_erf0_datt_erf_aux3*datt_erf_aux3_drho + datt_erf0_dityh_k_GGA*dityh_k_GGA_drho;

         double datt_erf0_dpi = datt_erf0_datt_erf_aux1*datt_erf_aux1_dpi + datt_erf0_datt_erf_aux2*datt_erf_aux2_dpi + datt_erf0_datt_erf_aux3*datt_erf_aux3_dpi + datt_erf0_dityh_k_GGA*dityh_k_GGA_dpi;

         datt_erf0_dsig = datt_erf0_datt_erf_aux1*datt_erf_aux1_dsig + datt_erf0_datt_erf_aux2*datt_erf_aux2_dsig + datt_erf0_datt_erf_aux3*datt_erf_aux3_dsig + datt_erf0_dityh_k_GGA*dityh_k_GGA_dsig;

         ityh_f = att_erf0/b88_denom;

         dityh_f_db88_denom = -att_erf0/pow(b88_denom, 2);

         dityh_f_datt_erf0 = 1.0/b88_denom;

         dityh_f_drho = datt_erf0_drho*dityh_f_datt_erf0 + db88_denom_drho*dityh_f_db88_denom;

         double dityh_f_dpi = datt_erf0_dpi*dityh_f_datt_erf0 + db88_denom_dpi*dityh_f_db88_denom;

         dityh_f_dsig = datt_erf0_dsig*dityh_f_datt_erf0 + db88_denom_dsig*dityh_f_db88_denom;

         f = ityh_f*pow(zeta_abs + 1, 2.0/3.0);

         double df_dzeta_abs = (2.0/3.0)*ityh_f/cbrt(zeta_abs + 1);

         double df_dityh_f = pow(zeta_abs + 1, 2.0/3.0);

         df_drho = 0;

         df_dpi = 0;

         df_dsig = 0;

         if (1 - zeta_abs > 1.0e-16)
         {
            double rho_s_2 = 0.5*density*(1 - zeta_abs);

            double drho_s_2_dzeta_abs = -0.5*density;

            double drho_s_2_drho = drho_s_2_dzeta_abs*dzeta_abs_drho - 0.5*zeta_abs + 0.5;

            double drho_s_2_dpi = drho_s_2_dzeta_abs*dzeta_abs_dpi;

            double xs_2 = cbrt(2)*xt/cbrt(1 - zeta_abs);

            double dxs_2_dxt = cbrt(2)/cbrt(1 - zeta_abs);

            double dxs_2_dzeta_abs = (1.0/3.0)*cbrt(2)*xt/pow(1 - zeta_abs, 4.0/3.0);

            double dxs_2_drho = dxs_2_dxt*dxt_drho + dxs_2_dzeta_abs*dzeta_abs_drho;

            double dxs_2_dsig = dxs_2_dxt*dxt_dsig;

            double dxs_2_dpi = dxs_2_dzeta_abs*dzeta_abs_dpi;

            double b88_denom_2 = b88_beta*b88_gamma*xs_2*asinh(xs_2) + 1;

            double db88_denom_2_dxs_2 = b88_beta*b88_gamma*xs_2/sqrt(pow(xs_2, 2) + 1) + b88_beta*b88_gamma*asinh(xs_2);

            double db88_denom_2_drho = db88_denom_2_dxs_2*dxs_2_drho;

            double db88_denom_2_dsig = db88_denom_2_dxs_2*dxs_2_dsig;

            double db88_denom_2_dpi = db88_denom_2_dxs_2*dxs_2_dpi;

            double b88_f_2 = b88_beta*pow(xs_2, 2)/(b88_denom_2*x_factor_c) + 1;

            double db88_f_2_dxs_2 = 2*b88_beta*xs_2/(b88_denom_2*x_factor_c);

            double db88_f_2_db88_denom_2 = -b88_beta*pow(xs_2, 2)/(pow(b88_denom_2, 2)*x_factor_c);

            double db88_f_2_drho = db88_denom_2_drho*db88_f_2_db88_denom_2 + db88_f_2_dxs_2*dxs_2_drho;

            double db88_f_2_dpi = db88_denom_2_dpi*db88_f_2_db88_denom_2 + db88_f_2_dxs_2*dxs_2_dpi;

            double db88_f_2_dsig = db88_denom_2_dsig*db88_f_2_db88_denom_2 + db88_f_2_dxs_2*dxs_2_dsig;

            double ityh_k_GGA_2 = (3.0/2.0)*M_SQRT2*sqrt(M_PI)*cbrt(rho_s_2)*sqrt(1/(b88_f_2*x_factor_c));

            double dityh_k_GGA_2_drho_s_2 = (1.0/2.0)*M_SQRT2*sqrt(M_PI)*sqrt(1/(b88_f_2*x_factor_c))/pow(rho_s_2, 2.0/3.0);

            double dityh_k_GGA_2_db88_f_2 = -3.0/4.0*M_SQRT2*sqrt(M_PI)*cbrt(rho_s_2)*sqrt(1/(b88_f_2*x_factor_c))/b88_f_2;

            double dityh_k_GGA_2_drho = db88_f_2_drho*dityh_k_GGA_2_db88_f_2 + dityh_k_GGA_2_drho_s_2*drho_s_2_drho;

            double dityh_k_GGA_2_dpi = db88_f_2_dpi*dityh_k_GGA_2_db88_f_2 + dityh_k_GGA_2_drho_s_2*drho_s_2_dpi;

            double dityh_k_GGA_2_dsig = db88_f_2_dsig*dityh_k_GGA_2_db88_f_2;

            double att_erf_aux1_2 = sqrt(M_PI)*erf(ityh_k_GGA_2/mu);

            double datt_erf_aux1_2_dityh_k_GGA_2 = 2*exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2))/mu;

            double datt_erf_aux1_2_drho = datt_erf_aux1_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

            double datt_erf_aux1_2_dsig = datt_erf_aux1_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

            double datt_erf_aux1_2_dpi = datt_erf_aux1_2_dityh_k_GGA_2*dityh_k_GGA_2_dpi;

            double datt_erf_aux2_2_dpi = 0;

            double datt_erf_aux2_2_dsig = 0;

            double att_erf_aux2_2 = 0;

            double datt_erf_aux2_2_dityh_k_GGA_2 = 0;

            double datt_erf_aux2_2_drho = 0;

            if ((1.0/2.0)*mu/ityh_k_GGA_2 < 5)
            {
               att_erf_aux2_2 = -1 + exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2));

               datt_erf_aux2_2_dityh_k_GGA_2 = -2*ityh_k_GGA_2*exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2))/pow(mu, 2);

               datt_erf_aux2_2_drho = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

               datt_erf_aux2_2_dsig = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

               datt_erf_aux2_2_dpi = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dpi;

            }
            else
            {
               att_erf_aux2_2 = -pow(ityh_k_GGA_2, 2)*(-2*pow(ityh_k_GGA_2, 2)/pow(mu, 2) + 1)/pow(mu, 2);

               datt_erf_aux2_2_dityh_k_GGA_2 = 4*pow(ityh_k_GGA_2, 3)/pow(mu, 4) - 2*ityh_k_GGA_2*(-2*pow(ityh_k_GGA_2, 2)/pow(mu, 2) + 1)/pow(mu, 2);

               datt_erf_aux2_2_drho = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

               datt_erf_aux2_2_dsig = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

               datt_erf_aux2_2_dpi = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dpi;

            }
            double att_erf_aux3_2 = (1.0/2.0)*att_erf_aux2_2*pow(mu, 2)/pow(ityh_k_GGA_2, 2) + 1.0/2.0;

            double datt_erf_aux3_2_dityh_k_GGA_2 = -att_erf_aux2_2*pow(mu, 2)/pow(ityh_k_GGA_2, 3);

            double datt_erf_aux3_2_datt_erf_aux2_2 = (1.0/2.0)*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

            double datt_erf_aux3_2_drho = datt_erf_aux2_2_drho*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

            double datt_erf_aux3_2_dsig = datt_erf_aux2_2_dsig*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

            double datt_erf_aux3_2_dpi = datt_erf_aux2_2_dpi*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_dityh_k_GGA_2*dityh_k_GGA_2_dpi;

            double att_erf0_2 = 1 - 4.0/3.0*mu*(att_erf_aux1_2 + mu*(att_erf_aux2_2 - att_erf_aux3_2)/ityh_k_GGA_2)/ityh_k_GGA_2;

            double datt_erf0_2_dityh_k_GGA_2 = (4.0/3.0)*mu*(att_erf_aux1_2 + mu*(att_erf_aux2_2 - att_erf_aux3_2)/ityh_k_GGA_2)/pow(ityh_k_GGA_2, 2) + (4.0/3.0)*pow(mu, 2)*(att_erf_aux2_2 - att_erf_aux3_2)/pow(ityh_k_GGA_2, 3);

            double datt_erf0_2_datt_erf_aux2_2 = -4.0/3.0*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

            double datt_erf0_2_datt_erf_aux3_2 = (4.0/3.0)*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

            double datt_erf0_2_datt_erf_aux1_2 = -4.0/3.0*mu/ityh_k_GGA_2;

            double datt_erf0_2_drho = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_drho + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_drho + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_drho + datt_erf0_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

            double datt_erf0_2_dpi = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_dpi + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_dpi + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_dpi + datt_erf0_2_dityh_k_GGA_2*dityh_k_GGA_2_dpi;

            double datt_erf0_2_dsig = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_dsig + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_dsig + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_dsig + datt_erf0_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

            double ityh_f_2 = att_erf0_2/b88_denom_2;

            double dityh_f_2_datt_erf0_2 = 1.0/b88_denom_2;

            double dityh_f_2_db88_denom_2 = -att_erf0_2/pow(b88_denom_2, 2);

            double dityh_f_2_drho = datt_erf0_2_drho*dityh_f_2_datt_erf0_2 + db88_denom_2_drho*dityh_f_2_db88_denom_2;

            double dityh_f_2_dpi = datt_erf0_2_dpi*dityh_f_2_datt_erf0_2 + db88_denom_2_dpi*dityh_f_2_db88_denom_2;

            double dityh_f_2_dsig = datt_erf0_2_dsig*dityh_f_2_datt_erf0_2 + db88_denom_2_dsig*dityh_f_2_db88_denom_2;

            f += ityh_f_2*pow(1 - zeta_abs, 2.0/3.0);

            double df_dityh_f_2 = pow(1 - zeta_abs, 2.0/3.0);

            df_dzeta_abs += - 2.0/3.0*ityh_f_2/cbrt(1 - zeta_abs);

            df_drho = df_dityh_f_2*dityh_f_2_drho;

            df_dpi = df_dityh_f_2*dityh_f_2_dpi;

            df_dsig = df_dityh_f_2*dityh_f_2_dsig;
         }

         df_drho += df_dityh_f*dityh_f_drho + df_dzeta_abs*dzeta_abs_drho;

         df_dpi += df_dityh_f*dityh_f_dpi + df_dzeta_abs*dzeta_abs_dpi;

         df_dsig += df_dityh_f*dityh_f_dsig;

      }
      else
      {
         // Lacks pi gradients
         rho_s = 0.5*density;

         drho_s_drho = 0.5;

         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         xs = cbrt(2)*xt;

         dxs_dxt = cbrt(2);

         dxs_drho = dxs_dxt*dxt_drho;

         dxs_dsig = dxs_dxt*dxt_dsig;

         b88_denom = b88_beta*b88_gamma*xs*asinh(xs) + 1;

         db88_denom_dxs = b88_beta*b88_gamma*xs/sqrt(pow(xs, 2) + 1) + b88_beta*b88_gamma*asinh(xs);

         db88_denom_drho = db88_denom_dxs*dxs_drho;

         db88_denom_dsig = db88_denom_dxs*dxs_dsig;

         b88_f = b88_beta*pow(xs, 2)/(b88_denom*x_factor_c) + 1;

         db88_f_db88_denom = -b88_beta*pow(xs, 2)/(pow(b88_denom, 2)*x_factor_c);

         db88_f_dxs = 2*b88_beta*xs/(b88_denom*x_factor_c);

         db88_f_drho = db88_denom_drho*db88_f_db88_denom + db88_f_dxs*dxs_drho;

         db88_f_dsig = db88_denom_dsig*db88_f_db88_denom + db88_f_dxs*dxs_dsig;

         ityh_k_GGA = (3.0/2.0)*M_SQRT2*sqrt(M_PI)*cbrt(rho_s)*sqrt(1/(b88_f*x_factor_c));

         dityh_k_GGA_db88_f = -3.0/4.0*M_SQRT2*sqrt(M_PI)*cbrt(rho_s)*sqrt(1/(b88_f*x_factor_c))/b88_f;

         dityh_k_GGA_drho = db88_f_drho*dityh_k_GGA_db88_f + (1.0/2.0)*M_SQRT2*sqrt(M_PI)*drho_s_drho*sqrt(1/(b88_f*x_factor_c))/pow(rho_s, 2.0/3.0);

         dityh_k_GGA_dsig = db88_f_dsig*dityh_k_GGA_db88_f;

         att_erf_aux1 = sqrt(M_PI)*erf(ityh_k_GGA/mu);

         datt_erf_aux1_dityh_k_GGA = 2*exp(-pow(ityh_k_GGA, 2)/pow(mu, 2))/mu;

         datt_erf_aux1_drho = datt_erf_aux1_dityh_k_GGA*dityh_k_GGA_drho;

         datt_erf_aux1_dsig = datt_erf_aux1_dityh_k_GGA*dityh_k_GGA_dsig;

         double datt_erf_aux2_dityh_k_GGA = 0;

         double att_erf_aux2 = 0;

         double datt_erf_aux2_dsig = 0;

         double datt_erf_aux2_drho = 0;

         if ((1.0/2.0)*mu/ityh_k_GGA < 5)
         {
            att_erf_aux2 = -1 + exp(-pow(ityh_k_GGA, 2)/pow(mu, 2));

            datt_erf_aux2_dityh_k_GGA = -2*ityh_k_GGA*exp(-pow(ityh_k_GGA, 2)/pow(mu, 2))/pow(mu, 2);

            datt_erf_aux2_drho = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_drho;

            datt_erf_aux2_dsig = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dsig;

         }
         else
         {
            att_erf_aux2 = -pow(ityh_k_GGA, 2)*(-2*pow(ityh_k_GGA, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_dityh_k_GGA = 4*pow(ityh_k_GGA, 3)/pow(mu, 4) - 2*ityh_k_GGA*(-2*pow(ityh_k_GGA, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_drho = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_drho;

            datt_erf_aux2_dsig = datt_erf_aux2_dityh_k_GGA*dityh_k_GGA_dsig;

         }
         att_erf_aux3 = (1.0/2.0)*att_erf_aux2*pow(mu, 2)/pow(ityh_k_GGA, 2) + 1.0/2.0;

         datt_erf_aux3_dityh_k_GGA = -att_erf_aux2*pow(mu, 2)/pow(ityh_k_GGA, 3);

         datt_erf_aux3_datt_erf_aux2 = (1.0/2.0)*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf_aux3_drho = datt_erf_aux2_drho*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_dityh_k_GGA*dityh_k_GGA_drho;

         datt_erf_aux3_dsig = datt_erf_aux2_dsig*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_dityh_k_GGA*dityh_k_GGA_dsig;

         att_erf0 = 1 - 4.0/3.0*mu*(att_erf_aux1 + mu*(att_erf_aux2 - att_erf_aux3)/ityh_k_GGA)/ityh_k_GGA;

         datt_erf0_dityh_k_GGA = (4.0/3.0)*mu*(att_erf_aux1 + mu*(att_erf_aux2 - att_erf_aux3)/ityh_k_GGA)/pow(ityh_k_GGA, 2) + (4.0/3.0)*pow(mu, 2)*(att_erf_aux2 - att_erf_aux3)/pow(ityh_k_GGA, 3);

         datt_erf0_datt_erf_aux2 = -4.0/3.0*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf0_datt_erf_aux3 = (4.0/3.0)*pow(mu, 2)/pow(ityh_k_GGA, 2);

         datt_erf0_datt_erf_aux1 = -4.0/3.0*mu/ityh_k_GGA;

         datt_erf0_drho = datt_erf0_datt_erf_aux1*datt_erf_aux1_drho + datt_erf0_datt_erf_aux2*datt_erf_aux2_drho + datt_erf0_datt_erf_aux3*datt_erf_aux3_drho + datt_erf0_dityh_k_GGA*dityh_k_GGA_drho;

         datt_erf0_dsig = datt_erf0_datt_erf_aux1*datt_erf_aux1_dsig + datt_erf0_datt_erf_aux2*datt_erf_aux2_dsig + datt_erf0_datt_erf_aux3*datt_erf_aux3_dsig + datt_erf0_dityh_k_GGA*dityh_k_GGA_dsig;

         ityh_f = att_erf0/b88_denom;

         dityh_f_db88_denom = -att_erf0/pow(b88_denom, 2);

         dityh_f_datt_erf0 = 1.0/b88_denom;

         dityh_f_drho = datt_erf0_drho*dityh_f_datt_erf0 + db88_denom_drho*dityh_f_db88_denom;

         dityh_f_dsig = datt_erf0_dsig*dityh_f_datt_erf0 + db88_denom_dsig*dityh_f_db88_denom;

         double f0 = 2*ityh_f;

         double df0_dityh_f = 2;

         double df0_drho = df0_dityh_f*dityh_f_drho;

         double df0_dsig = df0_dityh_f*dityh_f_dsig;

         double xs2 = pow(xs, 2);

         double dxs2_dxs = 2*xs;

         double dxs2_drho = dxs2_dxs*dxs_drho;

         double dxs2_dsig = dxs2_dxs*dxs_dsig;

         double sqrtx = sqrt(xs2 + 1);

         double dsqrtx_dxs2 = (1.0/2.0)/sqrt(xs2 + 1);

         double dsqrtx_drho = dsqrtx_dxs2*dxs2_drho;

         double dsqrtx_dsig = dsqrtx_dxs2*dxs2_dsig;

         double z2 = (1.0/9.0)*(-2*b88_beta*b88_gamma*xs2/sqrtx + b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 3) - 2*b88_denom + 2*pow(b88_beta*b88_gamma*xs2/sqrtx + b88_denom - 1, 2)/b88_denom)/pow(b88_denom, 2);

         double dz2_dsqrtx = (1.0/9.0)*(2*b88_beta*b88_gamma*xs2/pow(sqrtx, 2) - 3*b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 4) - 4*b88_beta*b88_gamma*xs2*(b88_beta*b88_gamma*xs2/sqrtx + b88_denom - 1)/(b88_denom*pow(sqrtx, 2)))/pow(b88_denom, 2);

         double dz2_db88_denom = (1.0/9.0)*(-2 + 2*(2*b88_beta*b88_gamma*xs2/sqrtx + 2*b88_denom - 2)/b88_denom - 2*pow(b88_beta*b88_gamma*xs2/sqrtx + b88_denom - 1, 2)/pow(b88_denom, 2))/pow(b88_denom, 2) - 2.0/9.0*(-2*b88_beta*b88_gamma*xs2/sqrtx + b88_beta*b88_gamma*pow(xs, 4)/pow(sqrtx, 3) - 2*b88_denom + 2*pow(b88_beta*b88_gamma*xs2/sqrtx + b88_denom - 1, 2)/b88_denom)/pow(b88_denom, 3);

         double dz2_dxs = (4.0/9.0)*b88_beta*b88_gamma*pow(xs, 3)/(pow(b88_denom, 2)*pow(sqrtx, 3));

         double dz2_dxs2 = (1.0/9.0)*(-2*b88_beta*b88_gamma/sqrtx + 4*b88_beta*b88_gamma*(b88_beta*b88_gamma*xs2/sqrtx + b88_denom - 1)/(b88_denom*sqrtx))/pow(b88_denom, 2);

         double dz2_drho = db88_denom_drho*dz2_db88_denom + dsqrtx_drho*dz2_dsqrtx + dxs2_drho*dz2_dxs2 + dxs_drho*dz2_dxs;

         double dz2_dsig = db88_denom_dsig*dz2_db88_denom + dsqrtx_dsig*dz2_dsqrtx + dxs2_dsig*dz2_dxs2 + dxs_dsig*dz2_dxs;

         double zeta2 = -2*pair_density/pow(density, 2);

         double dzeta2_drho = 4*pair_density/pow(density, 3);

         double dzeta2_dpi = -2/pow(density, 2);

         f = att_erf0*z2*zeta2 + f0;

         double df_df0 = 1;

         double df_dz2 = att_erf0*zeta2;

         double df_datt_erf0 = z2*zeta2;

         double df_dzeta2 = att_erf0*z2;

         df_drho = datt_erf0_drho*df_datt_erf0 + df0_drho*df_df0 + df_dz2*dz2_drho + df_dzeta2*dzeta2_drho;

         df_dsig = datt_erf0_dsig*df_datt_erf0 + df0_dsig*df_df0 + df_dz2*dz2_dsig;

         df_dpi = df_dzeta2*dzeta2_dpi;

      }
      double df2_drho = 0;

      double df2_dpi = 0;

      double f2 = 0;

      double df2_dsig = 0;

      double df2_df = 0;

      if (pair_density < 1.0000000000000001e-32)
      {
         f2 = f;

         df2_df = 1;

         df2_drho = df2_df*df_drho;

         df2_dsig = df2_df*df_dsig;

         df2_dpi = df2_df*df_dpi;

      }
      else
      {
         double rho_s_2 = 0.5*density;

         double drho_s_2_drho = 0.5;

         xt = pow(density, -1.3333333333333333)*sqrt(sig);

         dxt_drho = -1.3333333333333333*pow(density, -2.333333333333333)*sqrt(sig);

         dxt_dsig = (1.0/2.0)*pow(density, -1.3333333333333333)/sqrt(sig);

         double xs_2 = cbrt(2)*xt;

         double dxs_2_dxt = cbrt(2);

         double dxs_2_drho = dxs_2_dxt*dxt_drho;

         double dxs_2_dsig = dxs_2_dxt*dxt_dsig;

         double b88_denom_2 = b88_beta*b88_gamma*xs_2*asinh(xs_2) + 1;

         double db88_denom_2_dxs_2 = b88_beta*b88_gamma*xs_2/sqrt(pow(xs_2, 2) + 1) + b88_beta*b88_gamma*asinh(xs_2);

         double db88_denom_2_drho = db88_denom_2_dxs_2*dxs_2_drho;

         double db88_denom_2_dsig = db88_denom_2_dxs_2*dxs_2_dsig;

         double b88_f_2 = b88_beta*pow(xs_2, 2)/(b88_denom_2*x_factor_c) + 1;

         double db88_f_2_dxs_2 = 2*b88_beta*xs_2/(b88_denom_2*x_factor_c);

         double db88_f_2_db88_denom_2 = -b88_beta*pow(xs_2, 2)/(pow(b88_denom_2, 2)*x_factor_c);

         double db88_f_2_drho = db88_denom_2_drho*db88_f_2_db88_denom_2 + db88_f_2_dxs_2*dxs_2_drho;

         double db88_f_2_dsig = db88_denom_2_dsig*db88_f_2_db88_denom_2 + db88_f_2_dxs_2*dxs_2_dsig;

         double ityh_k_GGA_2 = (3.0/2.0)*M_SQRT2*sqrt(M_PI)*cbrt(rho_s_2)*sqrt(1/(b88_f_2*x_factor_c));

         double dityh_k_GGA_2_db88_f_2 = -3.0/4.0*M_SQRT2*sqrt(M_PI)*cbrt(rho_s_2)*sqrt(1/(b88_f_2*x_factor_c))/b88_f_2;

         double dityh_k_GGA_2_drho = db88_f_2_drho*dityh_k_GGA_2_db88_f_2 + (1.0/2.0)*M_SQRT2*sqrt(M_PI)*drho_s_2_drho*sqrt(1/(b88_f_2*x_factor_c))/pow(rho_s_2, 2.0/3.0);

         double dityh_k_GGA_2_dsig = db88_f_2_dsig*dityh_k_GGA_2_db88_f_2;

         double att_erf_aux1_2 = sqrt(M_PI)*erf(ityh_k_GGA_2/mu);

         double datt_erf_aux1_2_dityh_k_GGA_2 = 2*exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2))/mu;

         double datt_erf_aux1_2_drho = datt_erf_aux1_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

         double datt_erf_aux1_2_dsig = datt_erf_aux1_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

         double datt_erf_aux2_2_drho = 0;

         double datt_erf_aux2_2_dityh_k_GGA_2 = 0;

         double datt_erf_aux2_2_dsig = 0;

         double att_erf_aux2_2 = 0;

         if ((1.0/2.0)*mu/ityh_k_GGA_2 < 5)
         {
            att_erf_aux2_2 = -1 + exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2));

            datt_erf_aux2_2_dityh_k_GGA_2 = -2*ityh_k_GGA_2*exp(-pow(ityh_k_GGA_2, 2)/pow(mu, 2))/pow(mu, 2);

            datt_erf_aux2_2_drho = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

            datt_erf_aux2_2_dsig = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

         }
         else
         {
            att_erf_aux2_2 = -pow(ityh_k_GGA_2, 2)*(-2*pow(ityh_k_GGA_2, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_2_dityh_k_GGA_2 = 4*pow(ityh_k_GGA_2, 3)/pow(mu, 4) - 2*ityh_k_GGA_2*(-2*pow(ityh_k_GGA_2, 2)/pow(mu, 2) + 1)/pow(mu, 2);

            datt_erf_aux2_2_drho = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

            datt_erf_aux2_2_dsig = datt_erf_aux2_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

         }
         double att_erf_aux3_2 = (1.0/2.0)*att_erf_aux2_2*pow(mu, 2)/pow(ityh_k_GGA_2, 2) + 1.0/2.0;

         double datt_erf_aux3_2_dityh_k_GGA_2 = -att_erf_aux2_2*pow(mu, 2)/pow(ityh_k_GGA_2, 3);

         double datt_erf_aux3_2_datt_erf_aux2_2 = (1.0/2.0)*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

         double datt_erf_aux3_2_drho = datt_erf_aux2_2_drho*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

         double datt_erf_aux3_2_dsig = datt_erf_aux2_2_dsig*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

         double att_erf0_2 = 1 - 4.0/3.0*mu*(att_erf_aux1_2 + mu*(att_erf_aux2_2 - att_erf_aux3_2)/ityh_k_GGA_2)/ityh_k_GGA_2;

         double datt_erf0_2_dityh_k_GGA_2 = (4.0/3.0)*mu*(att_erf_aux1_2 + mu*(att_erf_aux2_2 - att_erf_aux3_2)/ityh_k_GGA_2)/pow(ityh_k_GGA_2, 2) + (4.0/3.0)*pow(mu, 2)*(att_erf_aux2_2 - att_erf_aux3_2)/pow(ityh_k_GGA_2, 3);

         double datt_erf0_2_datt_erf_aux2_2 = -4.0/3.0*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

         double datt_erf0_2_datt_erf_aux3_2 = (4.0/3.0)*pow(mu, 2)/pow(ityh_k_GGA_2, 2);

         double datt_erf0_2_datt_erf_aux1_2 = -4.0/3.0*mu/ityh_k_GGA_2;

         double datt_erf0_2_drho = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_drho + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_drho + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_drho + datt_erf0_2_dityh_k_GGA_2*dityh_k_GGA_2_drho;

         double datt_erf0_2_dsig = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_dsig + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_dsig + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_dsig + datt_erf0_2_dityh_k_GGA_2*dityh_k_GGA_2_dsig;

         double ityh_f_2 = att_erf0_2/b88_denom_2;

         double dityh_f_2_datt_erf0_2 = 1.0/b88_denom_2;

         double dityh_f_2_db88_denom_2 = -att_erf0_2/pow(b88_denom_2, 2);

         double dityh_f_2_drho = datt_erf0_2_drho*dityh_f_2_datt_erf0_2 + db88_denom_2_drho*dityh_f_2_db88_denom_2;

         double dityh_f_2_dsig = datt_erf0_2_dsig*dityh_f_2_datt_erf0_2 + db88_denom_2_dsig*dityh_f_2_db88_denom_2;

         f2 = -f + 4*ityh_f_2;

         double df2_dityh_f_2 = 4;

         df2_df = -1;

         df2_drho = df2_df*df_drho + df2_dityh_f_2*dityh_f_2_drho;

         df2_dpi = df2_df*df_dpi;

         df2_dsig = df2_df*df_dsig + df2_dityh_f_2*dityh_f_2_dsig;

      }
      double xt2 = pow(density, -2.6666666666666665)*sig;

      double dxt2_drho = -2.6666666666666665*pow(density, -3.6666666666666665)*sig;

      double dxt2_dsig = pow(density, -2.6666666666666665);

      exc[g] =  -1.0/2.0*cbrt(2)*b88_beta*cbrt(density)*f2*xt2;

      double dExc_df2 = -1.0/2.0*cbrt(2)*b88_beta*pow(density, 4.0/3.0)*xt2;

      double dExc_dxt2 = -1.0/2.0*cbrt(2)*b88_beta*pow(density, 4.0/3.0)*f2;

      vrho[2 * g + 0] =  -2.0/3.0*cbrt(2)*b88_beta*cbrt(density)*f2*xt2 + dExc_df2*df2_drho + dExc_dxt2*dxt2_drho;

      vrho[2 * g + 1] =  dExc_df2*df2_dpi;

      vsigma[3 * g + 0] =  dExc_df2*df2_dsig + dExc_dxt2*dxt2_dsig;

      //Currently, no explicit dependence
      vsigma[3 * g + 1] = 0.0;
      vsigma[3 * g + 2] = 0.0;

   }
}

}  // namespace pdftb88_erf