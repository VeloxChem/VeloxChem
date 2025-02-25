//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "PairDensityPBEX_erf.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftpbex_erf {  // pdftpbex_erf namespace

void
compute_exc_vxc(const int np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma, const double mu)
{

   // Subroutine generated by xc_write in MultiPsi, copyright M.G. Delcey, 2024

   double kappa = 0.80400000000000005;

   double pbe_b = 0.2195149727645171;

   double b0 = 7.0/81.0;

   double ax = 19.0;

   double a_cnst = (1.0/6.0)*pow(2, 2.0/3.0)*cbrt(3)*mu/cbrt(M_PI);

   double x2s = (1.0/12.0)*pow(6, 2.0/3.0)/pow(M_PI, 2.0/3.0);

   double rsfact = 0.90856029641606983*pow(M_PI, -0.33333333333333331);

   double fre = -3.0/8.0*cbrt(3)/cbrt(M_PI);

   for (int g = 0; g < np; g++)
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

      double rho83 = 0.15749013123685915*pow(density, 8.0/3.0);

      double drho83_drho = 0.41997368329829105*pow(density, 5.0/3.0);

      double xfact = 0.25*sig*pow(x2s, 2)/rho83;

      double dxfact_drho = -0.25*drho83_drho*sig*pow(x2s, 2)/pow(rho83, 2);

      double dxfact_dsig = 0.25*pow(x2s, 2)/rho83;

      double rs = rsfact/cbrt(density);

      double drs_drho = -1.0/3.0*rsfact/pow(density, 4.0/3.0);

      double a_fact = a_cnst*rs;

      double da_fact_drho = a_cnst*drs_drho;

      double mu_t = 0;

      double dc1_drho = 0;

      double c4 = 0;

      double dpbe_F_drho = 0;

      double dpbe_F_dsig = 0;

      double dkf_drho = 0;

      double dfpol_drho = 0;

      double c1 = 0;

      double att_erf_aux1 = 0;

      double datt_erf0_drho = 0;

      double dbmu_drho = 0;

      double bmu = 0;

      double erfmu = 0;

      double dmu_t4_drho = 0;

      double dc4_drho = 0;

      double att_erf0 = 0;

      double mu_t4 = 0;

      double dfpol_dsig = 0;

      double datt_erf_aux3_drho = 0;

      double datt_erf_aux1_drho = 0;

      double dc3_drho = 0;

      double kf = 0;

      double pbe_F = 0;

      double c3 = 0;

      double c2 = 0;

      double derfmu_drho = 0;

      double dmu_t2_drho = 0;

      double mu_t2 = 0;

      double dmu_t_drho = 0;

      double att_erf_aux3 = 0;

      double dc2_drho = 0;

      double fpol = 0;

      double dfpol_dpi = 0;

      if (fabs(pair_density) > 1.0000000000000001e-32)
      {
         double dzeta_abs_dpi = 0;

         double dzeta_abs_drho = 0;

         double zeta_abs = 0;

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
         double a = a_fact/cbrt(zeta_abs + 1);

         double da_dzeta_abs = -1.0/3.0*a_fact/pow(zeta_abs + 1, 4.0/3.0);

         double da_drho = da_dzeta_abs*dzeta_abs_drho + da_fact_drho/cbrt(zeta_abs + 1);

         double da_dpi = da_dzeta_abs*dzeta_abs_dpi;

         double rs_a = cbrt(2)*rs/cbrt(zeta_abs + 1);

         double drs_a_dzeta_abs = -1.0/3.0*cbrt(2)*rs/pow(zeta_abs + 1, 4.0/3.0);

         double drs_a_drho = drs_a_dzeta_abs*dzeta_abs_drho + cbrt(2)*drs_drho/cbrt(zeta_abs + 1);

         double drs_a_dpi = drs_a_dzeta_abs*dzeta_abs_dpi;

         double xsa2 = xfact/pow(zeta_abs + 1, 2.0/3.0);

         double dxsa2_dxfact = pow(zeta_abs + 1, -2.0/3.0);

         double dxsa2_dzeta_abs = -2.0/3.0*xfact/pow(zeta_abs + 1, 5.0/3.0);

         double dxsa2_drho = dxfact_drho*dxsa2_dxfact + dxsa2_dzeta_abs*dzeta_abs_drho;

         double dxsa2_dsig = dxfact_dsig*dxsa2_dxfact;

         double dxsa2_dpi = dxsa2_dzeta_abs*dzeta_abs_dpi;

         kf = cbrt(3)*pow(M_PI, 2.0/3.0)*cbrt(rs_a);

         double dkf_drs_a = (1.0/3.0)*cbrt(3)*pow(M_PI, 2.0/3.0)/pow(rs_a, 2.0/3.0);

         dkf_drho = dkf_drs_a*drs_a_drho;

         double dkf_dpi = dkf_drs_a*drs_a_dpi;

         mu_t = (1.0/2.0)*mu/kf;

         double dmu_t_dkf = -1.0/2.0*mu/pow(kf, 2);

         dmu_t_drho = dkf_drho*dmu_t_dkf;

         double dmu_t_dpi = dkf_dpi*dmu_t_dkf;

         mu_t2 = pow(mu_t, 2);

         double dmu_t2_dmu_t = 2*mu_t;

         dmu_t2_drho = dmu_t2_dmu_t*dmu_t_drho;

         double dmu_t2_dpi = dmu_t2_dmu_t*dmu_t_dpi;

         mu_t4 = pow(mu_t2, 2);

         double dmu_t4_dmu_t2 = 2*mu_t2;

         dmu_t4_drho = dmu_t2_drho*dmu_t4_dmu_t2;

         double dmu_t4_dpi = dmu_t2_dpi*dmu_t4_dmu_t2;

         erfmu = erf((1.0/2.0)/mu_t);

         double derfmu_dmu_t = -exp(-(1.0/4.0)/pow(mu_t, 2))/(sqrt(M_PI)*pow(mu_t, 2));

         derfmu_drho = derfmu_dmu_t*dmu_t_drho;

         double derfmu_dpi = derfmu_dmu_t*dmu_t_dpi;

         c1 = 22*mu_t2 + 144*mu_t4 + 1;

         double dc1_dmu_t2 = 22;

         double dc1_dmu_t4 = 144;

         dc1_drho = dc1_dmu_t2*dmu_t2_drho + dc1_dmu_t4*dmu_t4_drho;

         double dc1_dpi = dc1_dmu_t2*dmu_t2_dpi + dc1_dmu_t4*dmu_t4_dpi;

         c2 = 2*mu_t2*(72*mu_t2 - 7);

         double dc2_dmu_t2 = 288*mu_t2 - 14;

         dc2_drho = dc2_dmu_t2*dmu_t2_drho;

         double dc2_dpi = dc2_dmu_t2*dmu_t2_dpi;

         c3 = -864*mu_t4*(2*mu_t2 - 1);

         double dc3_dmu_t2 = -1728*mu_t4;

         double dc3_dmu_t4 = 864 - 1728*mu_t2;

         dc3_drho = dc3_dmu_t2*dmu_t2_drho + dc3_dmu_t4*dmu_t4_drho;

         double dc3_dpi = dc3_dmu_t2*dmu_t2_dpi + dc3_dmu_t4*dmu_t4_dpi;

         c4 = mu_t2*(8*sqrt(M_PI)*erfmu*mu_t - 24*mu_t2 + 32*mu_t4 - 3);

         double dc4_dmu_t4 = 32*mu_t2;

         double dc4_derfmu = 8*sqrt(M_PI)*mu_t*mu_t2;

         double dc4_dmu_t2 = 8*sqrt(M_PI)*erfmu*mu_t - 48*mu_t2 + 32*mu_t4 - 3;

         double dc4_dmu_t = 8*sqrt(M_PI)*erfmu*mu_t2;

         dc4_drho = dc4_derfmu*derfmu_drho + dc4_dmu_t*dmu_t_drho + dc4_dmu_t2*dmu_t2_drho + dc4_dmu_t4*dmu_t4_drho;

         double dc4_dpi = dc4_derfmu*derfmu_dpi + dc4_dmu_t*dmu_t_dpi + dc4_dmu_t2*dmu_t2_dpi + dc4_dmu_t4*dmu_t4_dpi;

         double dbt_dpi = 0;

         double bt = 0;

         double dbt_drho = 0;

         if (mu_t < 0.050000000000000003)
         {
            bt = (1.0/54.0)*c2/c4;

            double dbt_dc2 = (1.0/54.0)/c4;

            double dbt_dc4 = -1.0/54.0*c2/pow(c4, 2);

            dbt_drho = dbt_dc2*dc2_drho + dbt_dc4*dc4_drho;

            dbt_dpi = dbt_dc2*dc2_dpi + dbt_dc4*dc4_dpi;

         }
         else if (mu_t < 10000000000.0)
         {
            double expmu = exp((1.0/4.0)/mu_t2);

            double dexpmu_dmu_t2 = -1.0/4.0*exp((1.0/4.0)/mu_t2)/pow(mu_t2, 2);

            double dexpmu_drho = dexpmu_dmu_t2*dmu_t2_drho;

            double dexpmu_dpi = dexpmu_dmu_t2*dmu_t2_dpi;

            bt = (-c1 + c2*expmu)/(c3 + 54*c4*expmu);

            double dbt_dc3 = -(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2);

            double dbt_dc2 = expmu/(c3 + 54*c4*expmu);

            double dbt_dexpmu = c2/(c3 + 54*c4*expmu) - 54*c4*(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2);

            double dbt_dc1 = -1/(c3 + 54*c4*expmu);

            double dbt_dc4 = -54*expmu*(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2);

            dbt_drho = dbt_dc1*dc1_drho + dbt_dc2*dc2_drho + dbt_dc3*dc3_drho + dbt_dc4*dc4_drho + dbt_dexpmu*dexpmu_drho;

            dbt_dpi = dbt_dc1*dc1_dpi + dbt_dc2*dc2_dpi + dbt_dc3*dc3_dpi + dbt_dc4*dc4_dpi + dbt_dexpmu*dexpmu_dpi;

         }
         else
         {
            bt = -(1.0/17280.0)/mu_t4 + (1.0/72.0)/mu_t2 - (23.0/358400.0)/pow(mu_t, 6);

            double dbt_dmu_t2 = -(1.0/72.0)/pow(mu_t2, 2);

            double dbt_dmu_t4 = (1.0/17280.0)/pow(mu_t4, 2);

            double dbt_dmu_t = (69.0/179200.0)/pow(mu_t, 7);

            dbt_drho = dbt_dmu_t*dmu_t_drho + dbt_dmu_t2*dmu_t2_drho + dbt_dmu_t4*dmu_t4_drho;

            dbt_dpi = dbt_dmu_t*dmu_t_dpi + dbt_dmu_t2*dmu_t2_dpi + dbt_dmu_t4*dmu_t4_dpi;

         }
         bmu = bt*pbe_b*exp(-ax*pow(mu_t, 2))/b0;

         double dbmu_dbt = pbe_b*exp(-ax*pow(mu_t, 2))/b0;

         double dbmu_dmu_t = -2*ax*bt*mu_t*pbe_b*exp(-ax*pow(mu_t, 2))/b0;

         dbmu_drho = dbmu_dbt*dbt_drho + dbmu_dmu_t*dmu_t_drho;

         double dbmu_dpi = dbmu_dbt*dbt_dpi + dbmu_dmu_t*dmu_t_dpi;

         pbe_F = kappa - kappa/(bmu*xsa2/kappa + 1) + 1;

         double dpbe_F_dbmu = xsa2/pow(bmu*xsa2/kappa + 1, 2);

         double dpbe_F_dxsa2 = bmu/pow(bmu*xsa2/kappa + 1, 2);

         double dpbe_F_dpi = dbmu_dpi*dpbe_F_dbmu + dpbe_F_dxsa2*dxsa2_dpi;

         dpbe_F_drho = dbmu_drho*dpbe_F_dbmu + dpbe_F_dxsa2*dxsa2_drho;

         dpbe_F_dsig = dpbe_F_dxsa2*dxsa2_dsig;

         att_erf_aux1 = sqrt(M_PI)*erf((1.0/2.0)/a);

         double datt_erf_aux1_da = -exp(-(1.0/4.0)/pow(a, 2))/pow(a, 2);

         datt_erf_aux1_drho = da_drho*datt_erf_aux1_da;

         double datt_erf_aux1_dpi = da_dpi*datt_erf_aux1_da;

         double datt_erf_aux2_drho = 0;

         double att_erf_aux2 = 0;

         double datt_erf_aux2_da = 0;

         double datt_erf_aux2_dpi = 0;

         if (a < 5)
         {
            att_erf_aux2 = -1 + exp(-(1.0/4.0)/pow(a, 2));

            datt_erf_aux2_da = (1.0/2.0)*exp(-(1.0/4.0)/pow(a, 2))/pow(a, 3);

            datt_erf_aux2_drho = da_drho*datt_erf_aux2_da;

            datt_erf_aux2_dpi = da_dpi*datt_erf_aux2_da;

         }
         else
         {
            att_erf_aux2 = -1.0/4.0*(1 - (1.0/2.0)/pow(a, 2))/pow(a, 2);

            datt_erf_aux2_da = (1.0/2.0)*(1 - (1.0/2.0)/pow(a, 2))/pow(a, 3) - (1.0/4.0)/pow(a, 5);

            datt_erf_aux2_drho = da_drho*datt_erf_aux2_da;

            datt_erf_aux2_dpi = da_dpi*datt_erf_aux2_da;

         }
         att_erf_aux3 = 2*pow(a, 2)*att_erf_aux2 + 1.0/2.0;

         double datt_erf_aux3_da = 4*a*att_erf_aux2;

         double datt_erf_aux3_datt_erf_aux2 = 2*pow(a, 2);

         datt_erf_aux3_drho = da_drho*datt_erf_aux3_da + datt_erf_aux2_drho*datt_erf_aux3_datt_erf_aux2;

         double datt_erf_aux3_dpi = da_dpi*datt_erf_aux3_da + datt_erf_aux2_dpi*datt_erf_aux3_datt_erf_aux2;

         att_erf0 = -8.0/3.0*a*(2*a*(att_erf_aux2 - att_erf_aux3) + att_erf_aux1) + 1;

         double datt_erf0_datt_erf_aux3 = (16.0/3.0)*pow(a, 2);

         double datt_erf0_da = -16.0/3.0*a*(att_erf_aux2 - att_erf_aux3) - 8.0/3.0*a*(2*att_erf_aux2 - 2*att_erf_aux3) - 8.0/3.0*att_erf_aux1;

         double datt_erf0_datt_erf_aux2 = -16.0/3.0*pow(a, 2);

         double datt_erf0_datt_erf_aux1 = -8.0/3.0*a;

         datt_erf0_drho = da_drho*datt_erf0_da + datt_erf0_datt_erf_aux1*datt_erf_aux1_drho + datt_erf0_datt_erf_aux2*datt_erf_aux2_drho + datt_erf0_datt_erf_aux3*datt_erf_aux3_drho;

         double datt_erf0_dpi = da_dpi*datt_erf0_da + datt_erf0_datt_erf_aux1*datt_erf_aux1_dpi + datt_erf0_datt_erf_aux2*datt_erf_aux2_dpi + datt_erf0_datt_erf_aux3*datt_erf_aux3_dpi;

         double exa = att_erf0*pbe_F*pow(zeta_abs + 1, 4.0/3.0);

         double dexa_dpbe_F = att_erf0*pow(zeta_abs + 1, 4.0/3.0);

         double dexa_datt_erf0 = pbe_F*pow(zeta_abs + 1, 4.0/3.0);

         double dexa_dzeta_abs = (4.0/3.0)*att_erf0*pbe_F*cbrt(zeta_abs + 1);

         double dexa_dpi = datt_erf0_dpi*dexa_datt_erf0 + dexa_dpbe_F*dpbe_F_dpi + dexa_dzeta_abs*dzeta_abs_dpi;

         double dexa_drho = datt_erf0_drho*dexa_datt_erf0 + dexa_dpbe_F*dpbe_F_drho + dexa_dzeta_abs*dzeta_abs_drho;

         double dexa_dsig = dexa_dpbe_F*dpbe_F_dsig;

         fpol = exa;

         dfpol_dpi = dexa_dpi;

         dfpol_drho = dexa_drho;

         dfpol_dsig = dexa_dsig;

         if (1 - zeta_abs > 1.0e-16)
         {
            double b = a_fact/cbrt(1 - zeta_abs);

            double db_dzeta_abs = (1.0/3.0)*a_fact/pow(1 - zeta_abs, 4.0/3.0);

            double db_drho = da_fact_drho/cbrt(1 - zeta_abs) + db_dzeta_abs*dzeta_abs_drho;

            double db_dpi = db_dzeta_abs*dzeta_abs_dpi;

            double rs_b = cbrt(2)*rs/cbrt(1 - zeta_abs);

            double drs_b_dzeta_abs = (1.0/3.0)*cbrt(2)*rs/pow(1 - zeta_abs, 4.0/3.0);

            double drs_b_drho = drs_b_dzeta_abs*dzeta_abs_drho + cbrt(2)*drs_drho/cbrt(1 - zeta_abs);

            double drs_b_dpi = drs_b_dzeta_abs*dzeta_abs_dpi;

            double xsb2 = xfact/pow(1 - zeta_abs, 2.0/3.0);

            double dxsb2_dxfact = pow(1 - zeta_abs, -2.0/3.0);

            double dxsb2_dzeta_abs = (2.0/3.0)*xfact/pow(1 - zeta_abs, 5.0/3.0);

            double dxsb2_drho = dxfact_drho*dxsb2_dxfact + dxsb2_dzeta_abs*dzeta_abs_drho;

            double dxsb2_dsig = dxfact_dsig*dxsb2_dxfact;

            double dxsb2_dpi = dxsb2_dzeta_abs*dzeta_abs_dpi;

            double kf_2 = cbrt(3)*pow(M_PI, 2.0/3.0)*cbrt(rs_b);

            double dkf_2_drs_b = (1.0/3.0)*cbrt(3)*pow(M_PI, 2.0/3.0)/pow(rs_b, 2.0/3.0);

            double dkf_2_drho = dkf_2_drs_b*drs_b_drho;

            double dkf_2_dpi = dkf_2_drs_b*drs_b_dpi;

            double mu_t_2 = (1.0/2.0)*mu/kf_2;

            double dmu_t_2_dkf_2 = -1.0/2.0*mu/pow(kf_2, 2);

            double dmu_t_2_drho = dkf_2_drho*dmu_t_2_dkf_2;

            double dmu_t_2_dpi = dkf_2_dpi*dmu_t_2_dkf_2;

            double mu_t2_2 = pow(mu_t_2, 2);

            double dmu_t2_2_dmu_t_2 = 2*mu_t_2;

            double dmu_t2_2_drho = dmu_t2_2_dmu_t_2*dmu_t_2_drho;

            double dmu_t2_2_dpi = dmu_t2_2_dmu_t_2*dmu_t_2_dpi;

            double mu_t4_2 = pow(mu_t2_2, 2);

            double dmu_t4_2_dmu_t2_2 = 2*mu_t2_2;

            double dmu_t4_2_drho = dmu_t2_2_drho*dmu_t4_2_dmu_t2_2;

            double dmu_t4_2_dpi = dmu_t2_2_dpi*dmu_t4_2_dmu_t2_2;

            double erfmu_2 = erf((1.0/2.0)/mu_t_2);

            double derfmu_2_dmu_t_2 = -exp(-(1.0/4.0)/pow(mu_t_2, 2))/(sqrt(M_PI)*pow(mu_t_2, 2));

            double derfmu_2_drho = derfmu_2_dmu_t_2*dmu_t_2_drho;

            double derfmu_2_dpi = derfmu_2_dmu_t_2*dmu_t_2_dpi;

            double c1_2 = 22*mu_t2_2 + 144*mu_t4_2 + 1;

            double dc1_2_dmu_t2_2 = 22;

            double dc1_2_dmu_t4_2 = 144;

            double dc1_2_drho = dc1_2_dmu_t2_2*dmu_t2_2_drho + dc1_2_dmu_t4_2*dmu_t4_2_drho;

            double dc1_2_dpi = dc1_2_dmu_t2_2*dmu_t2_2_dpi + dc1_2_dmu_t4_2*dmu_t4_2_dpi;

            double c2_2 = 2*mu_t2_2*(72*mu_t2_2 - 7);

            double dc2_2_dmu_t2_2 = 288*mu_t2_2 - 14;

            double dc2_2_drho = dc2_2_dmu_t2_2*dmu_t2_2_drho;

            double dc2_2_dpi = dc2_2_dmu_t2_2*dmu_t2_2_dpi;

            double c3_2 = -864*mu_t4_2*(2*mu_t2_2 - 1);

            double dc3_2_dmu_t2_2 = -1728*mu_t4_2;

            double dc3_2_dmu_t4_2 = 864 - 1728*mu_t2_2;

            double dc3_2_drho = dc3_2_dmu_t2_2*dmu_t2_2_drho + dc3_2_dmu_t4_2*dmu_t4_2_drho;

            double dc3_2_dpi = dc3_2_dmu_t2_2*dmu_t2_2_dpi + dc3_2_dmu_t4_2*dmu_t4_2_dpi;

            double c4_2 = mu_t2_2*(8*sqrt(M_PI)*erfmu_2*mu_t_2 - 24*mu_t2_2 + 32*mu_t4_2 - 3);

            double dc4_2_derfmu_2 = 8*sqrt(M_PI)*mu_t2_2*mu_t_2;

            double dc4_2_dmu_t2_2 = 8*sqrt(M_PI)*erfmu_2*mu_t_2 - 48*mu_t2_2 + 32*mu_t4_2 - 3;

            double dc4_2_dmu_t_2 = 8*sqrt(M_PI)*erfmu_2*mu_t2_2;

            double dc4_2_dmu_t4_2 = 32*mu_t2_2;

            double dc4_2_drho = dc4_2_derfmu_2*derfmu_2_drho + dc4_2_dmu_t2_2*dmu_t2_2_drho + dc4_2_dmu_t4_2*dmu_t4_2_drho + dc4_2_dmu_t_2*dmu_t_2_drho;

            double dc4_2_dpi = dc4_2_derfmu_2*derfmu_2_dpi + dc4_2_dmu_t2_2*dmu_t2_2_dpi + dc4_2_dmu_t4_2*dmu_t4_2_dpi + dc4_2_dmu_t_2*dmu_t_2_dpi;

            double dbt_2_dpi = 0;

            double dbt_2_drho = 0;

            double bt_2 = 0;

            if (mu_t_2 < 0.050000000000000003)
            {
               bt_2 = (1.0/54.0)*c2_2/c4_2;

               double dbt_2_dc4_2 = -1.0/54.0*c2_2/pow(c4_2, 2);

               double dbt_2_dc2_2 = (1.0/54.0)/c4_2;

               dbt_2_drho = dbt_2_dc2_2*dc2_2_drho + dbt_2_dc4_2*dc4_2_drho;

               dbt_2_dpi = dbt_2_dc2_2*dc2_2_dpi + dbt_2_dc4_2*dc4_2_dpi;

            }
            else if (mu_t_2 < 10000000000.0)
            {
               double expmu = exp((1.0/4.0)/mu_t2_2);

               double dexpmu_dmu_t2_2 = -1.0/4.0*exp((1.0/4.0)/mu_t2_2)/pow(mu_t2_2, 2);

               double dexpmu_drho = dexpmu_dmu_t2_2*dmu_t2_2_drho;

               double dexpmu_dpi = dexpmu_dmu_t2_2*dmu_t2_2_dpi;

               bt_2 = (-c1_2 + c2_2*expmu)/(c3_2 + 54*c4_2*expmu);

               double dbt_2_dc3_2 = -(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2);

               double dbt_2_dc2_2 = expmu/(c3_2 + 54*c4_2*expmu);

               double dbt_2_dc1_2 = -1/(c3_2 + 54*c4_2*expmu);

               double dbt_2_dc4_2 = -54*expmu*(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2);

               double dbt_2_dexpmu = c2_2/(c3_2 + 54*c4_2*expmu) - 54*c4_2*(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2);

               dbt_2_drho = dbt_2_dc1_2*dc1_2_drho + dbt_2_dc2_2*dc2_2_drho + dbt_2_dc3_2*dc3_2_drho + dbt_2_dc4_2*dc4_2_drho + dbt_2_dexpmu*dexpmu_drho;

               dbt_2_dpi = dbt_2_dc1_2*dc1_2_dpi + dbt_2_dc2_2*dc2_2_dpi + dbt_2_dc3_2*dc3_2_dpi + dbt_2_dc4_2*dc4_2_dpi + dbt_2_dexpmu*dexpmu_dpi;

            }
            else
            {
               bt_2 = -(23.0/358400.0)/pow(mu_t_2, 6) - (1.0/17280.0)/mu_t4_2 + (1.0/72.0)/mu_t2_2;

               double dbt_2_dmu_t_2 = (69.0/179200.0)/pow(mu_t_2, 7);

               double dbt_2_dmu_t4_2 = (1.0/17280.0)/pow(mu_t4_2, 2);

               double dbt_2_dmu_t2_2 = -(1.0/72.0)/pow(mu_t2_2, 2);

               dbt_2_drho = dbt_2_dmu_t2_2*dmu_t2_2_drho + dbt_2_dmu_t4_2*dmu_t4_2_drho + dbt_2_dmu_t_2*dmu_t_2_drho;

               dbt_2_dpi = dbt_2_dmu_t2_2*dmu_t2_2_dpi + dbt_2_dmu_t4_2*dmu_t4_2_dpi + dbt_2_dmu_t_2*dmu_t_2_dpi;

            }
            double bmu_2 = bt_2*pbe_b*exp(-ax*pow(mu_t_2, 2))/b0;

            double dbmu_2_dmu_t_2 = -2*ax*bt_2*mu_t_2*pbe_b*exp(-ax*pow(mu_t_2, 2))/b0;

            double dbmu_2_dbt_2 = pbe_b*exp(-ax*pow(mu_t_2, 2))/b0;

            double dbmu_2_drho = dbmu_2_dbt_2*dbt_2_drho + dbmu_2_dmu_t_2*dmu_t_2_drho;

            double dbmu_2_dpi = dbmu_2_dbt_2*dbt_2_dpi + dbmu_2_dmu_t_2*dmu_t_2_dpi;

            double pbe_F_2 = kappa - kappa/(bmu_2*xsb2/kappa + 1) + 1;

            double dpbe_F_2_dbmu_2 = xsb2/pow(bmu_2*xsb2/kappa + 1, 2);

            double dpbe_F_2_dxsb2 = bmu_2/pow(bmu_2*xsb2/kappa + 1, 2);

            double dpbe_F_2_dpi = dbmu_2_dpi*dpbe_F_2_dbmu_2 + dpbe_F_2_dxsb2*dxsb2_dpi;

            double dpbe_F_2_drho = dbmu_2_drho*dpbe_F_2_dbmu_2 + dpbe_F_2_dxsb2*dxsb2_drho;

            double dpbe_F_2_dsig = dpbe_F_2_dxsb2*dxsb2_dsig;

            double att_erf_aux1_2 = sqrt(M_PI)*erf((1.0/2.0)/b);

            double datt_erf_aux1_2_db = -exp(-(1.0/4.0)/pow(b, 2))/pow(b, 2);

            double datt_erf_aux1_2_drho = datt_erf_aux1_2_db*db_drho;

            double datt_erf_aux1_2_dpi = datt_erf_aux1_2_db*db_dpi;

            double att_erf_aux2_2 = 0;

            double datt_erf_aux2_2_db = 0;

            double datt_erf_aux2_2_dpi = 0;

            double datt_erf_aux2_2_drho = 0;

            if (b < 5)
            {
               att_erf_aux2_2 = -1 + exp(-(1.0/4.0)/pow(b, 2));

               datt_erf_aux2_2_db = (1.0/2.0)*exp(-(1.0/4.0)/pow(b, 2))/pow(b, 3);

               datt_erf_aux2_2_drho = datt_erf_aux2_2_db*db_drho;

               datt_erf_aux2_2_dpi = datt_erf_aux2_2_db*db_dpi;

            }
            else
            {
               att_erf_aux2_2 = -1.0/4.0*(1 - (1.0/2.0)/pow(b, 2))/pow(b, 2);

               datt_erf_aux2_2_db = (1.0/2.0)*(1 - (1.0/2.0)/pow(b, 2))/pow(b, 3) - (1.0/4.0)/pow(b, 5);

               datt_erf_aux2_2_drho = datt_erf_aux2_2_db*db_drho;

               datt_erf_aux2_2_dpi = datt_erf_aux2_2_db*db_dpi;

            }
            double att_erf_aux3_2 = 2*att_erf_aux2_2*pow(b, 2) + 1.0/2.0;

            double datt_erf_aux3_2_db = 4*att_erf_aux2_2*b;

            double datt_erf_aux3_2_datt_erf_aux2_2 = 2*pow(b, 2);

            double datt_erf_aux3_2_drho = datt_erf_aux2_2_drho*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_db*db_drho;

            double datt_erf_aux3_2_dpi = datt_erf_aux2_2_dpi*datt_erf_aux3_2_datt_erf_aux2_2 + datt_erf_aux3_2_db*db_dpi;

            double att_erf0_2 = -8.0/3.0*b*(att_erf_aux1_2 + 2*b*(att_erf_aux2_2 - att_erf_aux3_2)) + 1;

            double datt_erf0_2_datt_erf_aux3_2 = (16.0/3.0)*pow(b, 2);

            double datt_erf0_2_db = -8.0/3.0*att_erf_aux1_2 - 16.0/3.0*b*(att_erf_aux2_2 - att_erf_aux3_2) - 8.0/3.0*b*(2*att_erf_aux2_2 - 2*att_erf_aux3_2);

            double datt_erf0_2_datt_erf_aux2_2 = -16.0/3.0*pow(b, 2);

            double datt_erf0_2_datt_erf_aux1_2 = -8.0/3.0*b;

            double datt_erf0_2_drho = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_drho + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_drho + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_drho + datt_erf0_2_db*db_drho;

            double datt_erf0_2_dpi = datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_dpi + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_dpi + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_dpi + datt_erf0_2_db*db_dpi;

            double exb = att_erf0_2*pbe_F_2*pow(1 - zeta_abs, 4.0/3.0);

            double dexb_datt_erf0_2 = pbe_F_2*pow(1 - zeta_abs, 4.0/3.0);

            double dexb_dzeta_abs = -4.0/3.0*att_erf0_2*pbe_F_2*cbrt(1 - zeta_abs);

            double dexb_dpbe_F_2 = att_erf0_2*pow(1 - zeta_abs, 4.0/3.0);

            double dexb_dpi = datt_erf0_2_dpi*dexb_datt_erf0_2 + dexb_dpbe_F_2*dpbe_F_2_dpi + dexb_dzeta_abs*dzeta_abs_dpi;

            double dexb_drho = datt_erf0_2_drho*dexb_datt_erf0_2 + dexb_dpbe_F_2*dpbe_F_2_drho + dexb_dzeta_abs*dzeta_abs_drho;

            double dexb_dsig = dexb_dpbe_F_2*dpbe_F_2_dsig;

            fpol += exb;

            dfpol_drho += dexb_drho;

            dfpol_dpi += dexb_dpi;

            dfpol_dsig += dexb_dsig;
         }
      }
      else
      {
         // Lacks pi gradient
         kf = 0.79370052598409979*cbrt(3)*pow(M_PI, 2.0/3.0)*cbrt(rs);

         dkf_drho = 0.26456684199469993*cbrt(3)*pow(M_PI, 2.0/3.0)*drs_drho/pow(rs, 2.0/3.0);

         mu_t = (1.0/2.0)*mu/kf;

         dmu_t_drho = -1.0/2.0*dkf_drho*mu/pow(kf, 2);

         mu_t2 = pow(mu_t, 2);

         dmu_t2_drho = 2*dmu_t_drho*mu_t;

         mu_t4 = pow(mu_t2, 2);

         dmu_t4_drho = 2*dmu_t2_drho*mu_t2;

         erfmu = erf((1.0/2.0)/mu_t);

         derfmu_drho = -dmu_t_drho*exp(-(1.0/4.0)/pow(mu_t, 2))/(sqrt(M_PI)*pow(mu_t, 2));

         c1 = 22*mu_t2 + 144*mu_t4 + 1;

         dc1_drho = 22*dmu_t2_drho + 144*dmu_t4_drho;

         c2 = 2*mu_t2*(72*mu_t2 - 7);

         dc2_drho = dmu_t2_drho*(288*mu_t2 - 14);

         c3 = -864*mu_t4*(2*mu_t2 - 1);

         dc3_drho = -1728*dmu_t2_drho*mu_t4 + dmu_t4_drho*(864 - 1728*mu_t2);

         c4 = mu_t2*(8*sqrt(M_PI)*erfmu*mu_t - 24*mu_t2 + 32*mu_t4 - 3);

         dc4_drho = 8*sqrt(M_PI)*derfmu_drho*mu_t*mu_t2 + dmu_t2_drho*(8*sqrt(M_PI)*erfmu*mu_t - 48*mu_t2 + 32*mu_t4 - 3) + 32*dmu_t4_drho*mu_t2 + 8*sqrt(M_PI)*dmu_t_drho*erfmu*mu_t2;

         double bt = 0;

         double dbt_drho = 0;

         if (mu_t < 0.050000000000000003)
         {
            bt = (1.0/54.0)*c2/c4;

            dbt_drho = -1.0/54.0*c2*dc4_drho/pow(c4, 2) + (1.0/54.0)*dc2_drho/c4;

         }
         else if (mu_t < 10000000000.0)
         {
            double expmu = exp((1.0/4.0)/mu_t2);

            double dexpmu_drho = -1.0/4.0*dmu_t2_drho*exp((1.0/4.0)/mu_t2)/pow(mu_t2, 2);

            bt = (-c1 + c2*expmu)/(c3 + 54*c4*expmu);

            dbt_drho = -dc1_drho/(c3 + 54*c4*expmu) + dc2_drho*expmu/(c3 + 54*c4*expmu) - dc3_drho*(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2) - 54*dc4_drho*expmu*(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2) + dexpmu_drho*(c2/(c3 + 54*c4*expmu) - 54*c4*(-c1 + c2*expmu)/pow(c3 + 54*c4*expmu, 2));

         }
         else
         {
            bt = -(1.0/17280.0)/mu_t4 + (1.0/72.0)/mu_t2 - (23.0/358400.0)/pow(mu_t, 6);

            dbt_drho = -1.0/72.0*dmu_t2_drho/pow(mu_t2, 2) + (1.0/17280.0)*dmu_t4_drho/pow(mu_t4, 2) + (69.0/179200.0)*dmu_t_drho/pow(mu_t, 7);

         }
         bmu = bt*pbe_b*exp(-ax*pow(mu_t, 2))/b0;

         dbmu_drho = -2*ax*bt*dmu_t_drho*mu_t*pbe_b*exp(-ax*pow(mu_t, 2))/b0 + dbt_drho*pbe_b*exp(-ax*pow(mu_t, 2))/b0;

         pbe_F = kappa - kappa/(bmu*xfact/kappa + 1) + 1;

         double dpbe_F_dxfact = bmu/pow(bmu*xfact/kappa + 1, 2);

         dpbe_F_drho = dbmu_drho*xfact/pow(bmu*xfact/kappa + 1, 2) + dpbe_F_dxfact*dxfact_drho;

         dpbe_F_dsig = dpbe_F_dxfact*dxfact_dsig;

         att_erf_aux1 = sqrt(M_PI)*erf((1.0/2.0)/a_fact);

         datt_erf_aux1_drho = -da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 2);

         double datt_erf_aux2_drho = 0;

         double att_erf_aux2 = 0;

         if (a_fact < 5)
         {
            att_erf_aux2 = -1 + exp(-(1.0/4.0)/pow(a_fact, 2));

            datt_erf_aux2_drho = (1.0/2.0)*da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 3);

         }
         else
         {
            att_erf_aux2 = -1.0/4.0*(1 - (1.0/2.0)/pow(a_fact, 2))/pow(a_fact, 2);

            datt_erf_aux2_drho = da_fact_drho*((1.0/2.0)*(1 - (1.0/2.0)/pow(a_fact, 2))/pow(a_fact, 3) - (1.0/4.0)/pow(a_fact, 5));

         }
         att_erf_aux3 = 2*pow(a_fact, 2)*att_erf_aux2 + 1.0/2.0;

         datt_erf_aux3_drho = 2*pow(a_fact, 2)*datt_erf_aux2_drho + 4*a_fact*att_erf_aux2*da_fact_drho;

         att_erf0 = -8.0/3.0*a_fact*(2*a_fact*(att_erf_aux2 - att_erf_aux3) + att_erf_aux1) + 1;

         datt_erf0_drho = -16.0/3.0*pow(a_fact, 2)*datt_erf_aux2_drho + (16.0/3.0)*pow(a_fact, 2)*datt_erf_aux3_drho - 8.0/3.0*a_fact*datt_erf_aux1_drho + da_fact_drho*(-16.0/3.0*a_fact*(att_erf_aux2 - att_erf_aux3) - 8.0/3.0*a_fact*(2*att_erf_aux2 - 2*att_erf_aux3) - 8.0/3.0*att_erf_aux1);

         fpol = 2.0*att_erf0*pbe_F;

         double dfpol_dpbe_F = 2.0*att_erf0;

         dfpol_drho = 2.0*datt_erf0_drho*pbe_F + dfpol_dpbe_F*dpbe_F_drho;

         dfpol_dsig = dfpol_dpbe_F*dpbe_F_dsig;

      }
      double df_drho = 0;

      double f = 0;

      double df_dsig = 0;

      double df_dfpol = 0;

      double df_dpi = 0;

      if (pair_density < 1.0000000000000001e-32)
      {
         f = fpol;

         df_drho = dfpol_drho;

         df_dsig = dfpol_dsig;

         df_dpi = dfpol_dpi;

      }
      else
      {
         double kf_2 = 0.79370052598409979*cbrt(3)*pow(M_PI, 2.0/3.0)*cbrt(rs);

         double dkf_2_drho = 0.26456684199469993*cbrt(3)*pow(M_PI, 2.0/3.0)*drs_drho/pow(rs, 2.0/3.0);

         double mu_t_2 = (1.0/2.0)*mu/kf_2;

         double dmu_t_2_drho = -1.0/2.0*dkf_2_drho*mu/pow(kf_2, 2);

         double mu_t2_2 = pow(mu_t_2, 2);

         double dmu_t2_2_drho = 2*dmu_t_2_drho*mu_t_2;

         double mu_t4_2 = pow(mu_t2_2, 2);

         double dmu_t4_2_drho = 2*dmu_t2_2_drho*mu_t2_2;

         double erfmu_2 = erf((1.0/2.0)/mu_t_2);

         double derfmu_2_drho = -dmu_t_2_drho*exp(-(1.0/4.0)/pow(mu_t_2, 2))/(sqrt(M_PI)*pow(mu_t_2, 2));

         double c1_2 = 22*mu_t2_2 + 144*mu_t4_2 + 1;

         double dc1_2_drho = 22*dmu_t2_2_drho + 144*dmu_t4_2_drho;

         double c2_2 = 2*mu_t2_2*(72*mu_t2_2 - 7);

         double dc2_2_drho = dmu_t2_2_drho*(288*mu_t2_2 - 14);

         double c3_2 = -864*mu_t4_2*(2*mu_t2_2 - 1);

         double dc3_2_drho = -1728*dmu_t2_2_drho*mu_t4_2 + dmu_t4_2_drho*(864 - 1728*mu_t2_2);

         double c4_2 = mu_t2_2*(8*sqrt(M_PI)*erfmu_2*mu_t_2 - 24*mu_t2_2 + 32*mu_t4_2 - 3);

         double dc4_2_drho = 8*sqrt(M_PI)*derfmu_2_drho*mu_t2_2*mu_t_2 + dmu_t2_2_drho*(8*sqrt(M_PI)*erfmu_2*mu_t_2 - 48*mu_t2_2 + 32*mu_t4_2 - 3) + 32*dmu_t4_2_drho*mu_t2_2 + 8*sqrt(M_PI)*dmu_t_2_drho*erfmu_2*mu_t2_2;

         double dbt_2_drho = 0;

         double bt_2 = 0;

         if (mu_t_2 < 0.050000000000000003)
         {
            bt_2 = (1.0/54.0)*c2_2/c4_2;

            dbt_2_drho = -1.0/54.0*c2_2*dc4_2_drho/pow(c4_2, 2) + (1.0/54.0)*dc2_2_drho/c4_2;

         }
         else if (mu_t_2 < 10000000000.0)
         {
            double expmu = exp((1.0/4.0)/mu_t2_2);

            double dexpmu_drho = -1.0/4.0*dmu_t2_2_drho*exp((1.0/4.0)/mu_t2_2)/pow(mu_t2_2, 2);

            bt_2 = (-c1_2 + c2_2*expmu)/(c3_2 + 54*c4_2*expmu);

            dbt_2_drho = -dc1_2_drho/(c3_2 + 54*c4_2*expmu) + dc2_2_drho*expmu/(c3_2 + 54*c4_2*expmu) - dc3_2_drho*(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2) - 54*dc4_2_drho*expmu*(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2) + dexpmu_drho*(c2_2/(c3_2 + 54*c4_2*expmu) - 54*c4_2*(-c1_2 + c2_2*expmu)/pow(c3_2 + 54*c4_2*expmu, 2));

         }
         else
         {
            bt_2 = -(23.0/358400.0)/pow(mu_t_2, 6) - (1.0/17280.0)/mu_t4_2 + (1.0/72.0)/mu_t2_2;

            dbt_2_drho = -1.0/72.0*dmu_t2_2_drho/pow(mu_t2_2, 2) + (1.0/17280.0)*dmu_t4_2_drho/pow(mu_t4_2, 2) + (69.0/179200.0)*dmu_t_2_drho/pow(mu_t_2, 7);

         }
         double bmu_2 = bt_2*pbe_b*exp(-ax*pow(mu_t_2, 2))/b0;

         double dbmu_2_drho = -2*ax*bt_2*dmu_t_2_drho*mu_t_2*pbe_b*exp(-ax*pow(mu_t_2, 2))/b0 + dbt_2_drho*pbe_b*exp(-ax*pow(mu_t_2, 2))/b0;

         double pbe_F_2 = kappa - kappa/(bmu_2*xfact/kappa + 1) + 1;

         double dpbe_F_2_dxfact = bmu_2/pow(bmu_2*xfact/kappa + 1, 2);

         double dpbe_F_2_drho = dbmu_2_drho*xfact/pow(bmu_2*xfact/kappa + 1, 2) + dpbe_F_2_dxfact*dxfact_drho;

         double dpbe_F_2_dsig = dpbe_F_2_dxfact*dxfact_dsig;

         double att_erf_aux1_2 = sqrt(M_PI)*erf((1.0/2.0)/a_fact);

         double datt_erf_aux1_2_drho = -da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 2);

         double att_erf_aux2_2 = 0;

         double datt_erf_aux2_2_drho = 0;

         if (a_fact < 5)
         {
            att_erf_aux2_2 = -1 + exp(-(1.0/4.0)/pow(a_fact, 2));

            datt_erf_aux2_2_drho = (1.0/2.0)*da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 3);

         }
         else
         {
            att_erf_aux2_2 = -1.0/4.0*(1 - (1.0/2.0)/pow(a_fact, 2))/pow(a_fact, 2);

            datt_erf_aux2_2_drho = da_fact_drho*((1.0/2.0)*(1 - (1.0/2.0)/pow(a_fact, 2))/pow(a_fact, 3) - (1.0/4.0)/pow(a_fact, 5));

         }
         double att_erf_aux3_2 = 2*pow(a_fact, 2)*att_erf_aux2_2 + 1.0/2.0;

         double datt_erf_aux3_2_drho = 2*pow(a_fact, 2)*datt_erf_aux2_2_drho + 4*a_fact*att_erf_aux2_2*da_fact_drho;

         double att_erf0_2 = -8.0/3.0*a_fact*(2*a_fact*(att_erf_aux2_2 - att_erf_aux3_2) + att_erf_aux1_2) + 1;

         double datt_erf0_2_drho = -16.0/3.0*pow(a_fact, 2)*datt_erf_aux2_2_drho + (16.0/3.0)*pow(a_fact, 2)*datt_erf_aux3_2_drho - 8.0/3.0*a_fact*datt_erf_aux1_2_drho + da_fact_drho*(-16.0/3.0*a_fact*(att_erf_aux2_2 - att_erf_aux3_2) - 8.0/3.0*a_fact*(2*att_erf_aux2_2 - 2*att_erf_aux3_2) - 8.0/3.0*att_erf_aux1_2);

         f = 4.0*att_erf0_2*pbe_F_2 - fpol;

         df_dfpol = -1;

         double df_dpbe_F_2 = 4.0*att_erf0_2;

         df_drho = 4.0*datt_erf0_2_drho*pbe_F_2 + df_dfpol*dfpol_drho + df_dpbe_F_2*dpbe_F_2_drho;

         df_dpi = df_dfpol*dfpol_dpi;

         df_dsig = df_dfpol*dfpol_dsig + df_dpbe_F_2*dpbe_F_2_dsig;

      }
      exc[g] =  cbrt(density)*f*fre;

      double dExc_df = pow(density, 4.0/3.0)*fre;

      vrho[2 * g + 0] =  dExc_df*df_drho + (4.0/3.0)*cbrt(density)*f*fre;

      vrho[2 * g + 1] =  dExc_df*df_dpi;

      vsigma[3 * g + 0] =  dExc_df*df_dsig;

      //Currently, no explicit dependence
      vsigma[3 * g + 1] = 0.0;
      vsigma[3 * g + 2] = 0.0;

   }
}

}  // namespace pdftpbex_erf
