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

#include "PairDensitySlater_erf.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace pdftslater_erf {  // pdftslater_erf namespace

void
compute_exc_vxc(const int np, const double* rho, double* exc, double* vrho, const double mu)
{

   // Subroutine generated by xc_write in MultiPsi, copyright M.G. Delcey, 2024

   double a_cnst = (1.0/6.0)*pow(2, 2.0/3.0)*cbrt(3)*mu/cbrt(M_PI);

   double rs_factor = (1.0/2.0)*cbrt(6)/cbrt(M_PI);

   double x_factor_c = (3.0/4.0)*cbrt(6)/cbrt(M_PI);

   double lda_x_ax = -1.0/4.0*pow(2, 2.0/3.0)*rs_factor*x_factor_c;

   double rsfact = 0.90856029641606983*pow(M_PI, -0.33333333333333331);

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

      double pair_density = rho[2 * g + 1];

      double rs = rsfact/cbrt(density);

      double drs_drho = -1.0/3.0*rsfact/pow(density, 4.0/3.0);

      double a_fact = a_cnst*rs;

      double da_fact_drho = a_cnst*drs_drho;

      double datt_erf_aux3_drho = 0;

      double datt_erf0_drho = 0;

      double att_erf_aux1 = 0;

      double datt_erf_aux1_drho = 0;

      double dfpol_dpi = 0;

      double att_erf0 = 0;

      double dfpol_drho = 0;

      double fpol = 0;

      double att_erf_aux3 = 0;

      if (fabs(pair_density) > 9.9999999999999998e-17)
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
         double a = a_fact/cbrt(zeta_abs + 1);

         double da_dzeta_abs = -1.0/3.0*a_fact/pow(zeta_abs + 1, 4.0/3.0);

         double da_drho = da_dzeta_abs*dzeta_abs_drho + da_fact_drho/cbrt(zeta_abs + 1);

         double da_dpi = da_dzeta_abs*dzeta_abs_dpi;

         double att_erf_aux1_2 = sqrt(M_PI)*erf((1.0/2.0)/a);

         double datt_erf_aux1_2_da = -exp(-(1.0/4.0)/pow(a, 2))/pow(a, 2);

         double datt_erf_aux1_2_drho = da_drho*datt_erf_aux1_2_da;

         double datt_erf_aux1_2_dpi = da_dpi*datt_erf_aux1_2_da;

         double datt_erf_aux2_2_da = 0;

         double datt_erf_aux2_2_drho = 0;

         double att_erf_aux2_2 = 0;

         double datt_erf_aux2_2_dpi = 0;

         if (a < 5)
         {
            att_erf_aux2_2 = -1 + exp(-(1.0/4.0)/pow(a, 2));

            datt_erf_aux2_2_da = (1.0/2.0)*exp(-(1.0/4.0)/pow(a, 2))/pow(a, 3);

            datt_erf_aux2_2_drho = da_drho*datt_erf_aux2_2_da;

            datt_erf_aux2_2_dpi = da_dpi*datt_erf_aux2_2_da;

         }
         else
         {
            att_erf_aux2_2 = -1.0/4.0*(1 - (1.0/2.0)/pow(a, 2))/pow(a, 2);

            datt_erf_aux2_2_da = (1.0/2.0)*(1 - (1.0/2.0)/pow(a, 2))/pow(a, 3) - (1.0/4.0)/pow(a, 5);

            datt_erf_aux2_2_drho = da_drho*datt_erf_aux2_2_da;

            datt_erf_aux2_2_dpi = da_dpi*datt_erf_aux2_2_da;

         }
         double att_erf_aux3_2 = 2*pow(a, 2)*att_erf_aux2_2 + 1.0/2.0;

         double datt_erf_aux3_2_datt_erf_aux2_2 = 2*pow(a, 2);

         double datt_erf_aux3_2_da = 4*a*att_erf_aux2_2;

         double datt_erf_aux3_2_drho = da_drho*datt_erf_aux3_2_da + datt_erf_aux2_2_drho*datt_erf_aux3_2_datt_erf_aux2_2;

         double datt_erf_aux3_2_dpi = da_dpi*datt_erf_aux3_2_da + datt_erf_aux2_2_dpi*datt_erf_aux3_2_datt_erf_aux2_2;

         double att_erf0_2 = -8.0/3.0*a*(2*a*(att_erf_aux2_2 - att_erf_aux3_2) + att_erf_aux1_2) + 1;

         double datt_erf0_2_datt_erf_aux2_2 = -16.0/3.0*pow(a, 2);

         double datt_erf0_2_da = -16.0/3.0*a*(att_erf_aux2_2 - att_erf_aux3_2) - 8.0/3.0*a*(2*att_erf_aux2_2 - 2*att_erf_aux3_2) - 8.0/3.0*att_erf_aux1_2;

         double datt_erf0_2_datt_erf_aux1_2 = -8.0/3.0*a;

         double datt_erf0_2_datt_erf_aux3_2 = (16.0/3.0)*pow(a, 2);

         double datt_erf0_2_drho = da_drho*datt_erf0_2_da + datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_drho + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_drho + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_drho;

         double datt_erf0_2_dpi = da_dpi*datt_erf0_2_da + datt_erf0_2_datt_erf_aux1_2*datt_erf_aux1_2_dpi + datt_erf0_2_datt_erf_aux2_2*datt_erf_aux2_2_dpi + datt_erf0_2_datt_erf_aux3_2*datt_erf_aux3_2_dpi;

         fpol = att_erf0_2*pow(zeta_abs + 1, 4.0/3.0);

         double dfpol_datt_erf0_2 = pow(zeta_abs + 1, 4.0/3.0);

         double dfpol_dzeta_abs = (4.0/3.0)*att_erf0_2*cbrt(zeta_abs + 1);

         dfpol_drho = datt_erf0_2_drho*dfpol_datt_erf0_2;

         dfpol_dpi = datt_erf0_2_dpi*dfpol_datt_erf0_2;

         if ( 1-zeta_abs > 1.0e-16)
         {
             double b = a_fact/cbrt(1 - zeta_abs);

             double db_dzeta_abs = (1.0/3.0)*a_fact/pow(1 - zeta_abs, 4.0/3.0);

             double db_drho = da_fact_drho/cbrt(1 - zeta_abs) + db_dzeta_abs*dzeta_abs_drho;

             double db_dpi = db_dzeta_abs*dzeta_abs_dpi;

             att_erf_aux1 = sqrt(M_PI)*erf((1.0/2.0)/b);

             double datt_erf_aux1_db = -exp(-(1.0/4.0)/pow(b, 2))/pow(b, 2);

             datt_erf_aux1_drho = datt_erf_aux1_db*db_drho;

             double datt_erf_aux1_dpi = datt_erf_aux1_db*db_dpi;

             double att_erf_aux2 = 0;

             double datt_erf_aux2_dpi = 0;

             double datt_erf_aux2_db = 0;

             double datt_erf_aux2_drho = 0;

             if (b < 5)
             {
                att_erf_aux2 = -1 + exp(-(1.0/4.0)/pow(b, 2));

                datt_erf_aux2_db = (1.0/2.0)*exp(-(1.0/4.0)/pow(b, 2))/pow(b, 3);

                datt_erf_aux2_drho = datt_erf_aux2_db*db_drho;

                datt_erf_aux2_dpi = datt_erf_aux2_db*db_dpi;

             }
             else
             {
                att_erf_aux2 = -1.0/4.0*(1 - (1.0/2.0)/pow(b, 2))/pow(b, 2);

                datt_erf_aux2_db = (1.0/2.0)*(1 - (1.0/2.0)/pow(b, 2))/pow(b, 3) - (1.0/4.0)/pow(b, 5);

                datt_erf_aux2_drho = datt_erf_aux2_db*db_drho;

                datt_erf_aux2_dpi = datt_erf_aux2_db*db_dpi;

             }
             att_erf_aux3 = 2*att_erf_aux2*pow(b, 2) + 1.0/2.0;

             double datt_erf_aux3_db = 4*att_erf_aux2*b;

             double datt_erf_aux3_datt_erf_aux2 = 2*pow(b, 2);

             datt_erf_aux3_drho = datt_erf_aux2_drho*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_db*db_drho;

             double datt_erf_aux3_dpi = datt_erf_aux2_dpi*datt_erf_aux3_datt_erf_aux2 + datt_erf_aux3_db*db_dpi;

             att_erf0 = -8.0/3.0*b*(att_erf_aux1 + 2*b*(att_erf_aux2 - att_erf_aux3)) + 1;

             double datt_erf0_db = -8.0/3.0*att_erf_aux1 - 16.0/3.0*b*(att_erf_aux2 - att_erf_aux3) - 8.0/3.0*b*(2*att_erf_aux2 - 2*att_erf_aux3);

             double datt_erf0_datt_erf_aux2 = -16.0/3.0*pow(b, 2);

             double datt_erf0_datt_erf_aux3 = (16.0/3.0)*pow(b, 2);

             double datt_erf0_datt_erf_aux1 = -8.0/3.0*b;

             datt_erf0_drho = datt_erf0_datt_erf_aux1*datt_erf_aux1_drho + datt_erf0_datt_erf_aux2*datt_erf_aux2_drho + datt_erf0_datt_erf_aux3*datt_erf_aux3_drho + datt_erf0_db*db_drho;

             double datt_erf0_dpi = datt_erf0_datt_erf_aux1*datt_erf_aux1_dpi + datt_erf0_datt_erf_aux2*datt_erf_aux2_dpi + datt_erf0_datt_erf_aux3*datt_erf_aux3_dpi + datt_erf0_db*db_dpi;

             fpol += att_erf0*pow(1 - zeta_abs, 4.0/3.0);

             dfpol_dzeta_abs += -4.0/3.0*att_erf0*cbrt(1 - zeta_abs);

             double dfpol_datt_erf0 = pow(1 - zeta_abs, 4.0/3.0);

             dfpol_drho += datt_erf0_drho*dfpol_datt_erf0;

             dfpol_dpi += datt_erf0_dpi*dfpol_datt_erf0;
         }

         dfpol_drho += dfpol_dzeta_abs*dzeta_abs_drho;

         dfpol_dpi += dfpol_dzeta_abs*dzeta_abs_dpi;
      }
      else
      {
         double expa = exp(-(1.0/4.0)/pow(a_fact, 2));

         double dexpa_drho = (1.0/2.0)*da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 3);

         att_erf_aux1 = sqrt(M_PI)*erf((1.0/2.0)/a_fact);

         datt_erf_aux1_drho = -da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 2);

         double att_erf_aux2 = 0;

         double datt_erf_aux2_drho = 0;

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

         double zeta2 = -2*pair_density/pow(density, 2);

         double dzeta2_drho = 4*pair_density/pow(density, 3);

         double dzeta2_dpi = -2/pow(density, 2);

         fpol = 2.0*att_erf0 + 0.44444444444444442*zeta2*(4.0*pow(a_fact, 2)*(expa - 1) + 1.0);

         double dfpol_dzeta2 = 1.7777777777777777*pow(a_fact, 2)*(expa - 1) + 0.44444444444444442;

         dfpol_drho = 1.7777777777777777*pow(a_fact, 2)*dexpa_drho*zeta2 + 3.5555555555555554*a_fact*da_fact_drho*zeta2*(expa - 1) + 2.0*datt_erf0_drho + dfpol_dzeta2*dzeta2_drho;

         dfpol_dpi = dfpol_dzeta2*dzeta2_dpi;

      }
      double f = 0;

      double df_dpi = 0;

      double df_drho = 0;

      double df_dfpol = 0;

      if (pair_density < 9.9999999999999998e-17)
      {
         f = fpol;

         df_dfpol = 1;

         df_drho = df_dfpol*dfpol_drho;

         df_dpi = df_dfpol*dfpol_dpi;

      }
      else
      {
         double att_erf_aux1_2 = sqrt(M_PI)*erf((1.0/2.0)/a_fact);

         double datt_erf_aux1_2_drho = -da_fact_drho*exp(-(1.0/4.0)/pow(a_fact, 2))/pow(a_fact, 2);

         double datt_erf_aux2_2_drho = 0;

         double att_erf_aux2_2 = 0;

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

         f = 4.0*att_erf0_2 - fpol;

         df_dfpol = -1;

         df_drho = 4.0*datt_erf0_2_drho + df_dfpol*dfpol_drho;

         df_dpi = df_dfpol*dfpol_dpi;

      }
      exc[g] =  f*lda_x_ax/rs;

      double dExc_df = density*lda_x_ax/rs;

      vrho[2 * g + 0] =  dExc_df*df_drho - density*drs_drho*f*lda_x_ax/pow(rs, 2) + f*lda_x_ax/rs;

      vrho[2 * g + 1] =  dExc_df*df_dpi;

   }
}

}  // namespace pdftslater_erf


