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

#include "PkzbFunctional.hpp"

#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace vxcfuncs {  // vxcfuncs namespace

    CXCFunctional
    setPkzbFunctional()
    {
        return CXCFunctional({"PKZB"}, xcfun::mgga, 0.0, {setPrimitivePkzbFunctional()}, {1.0});
    }

    CPrimitiveFunctional
    setPrimitivePkzbFunctional()
    {
        return CPrimitiveFunctional({"PKZB"},
                                    xcfun::mgga,
                                    &vxcfuncs::PkzbFuncGradientAB,
                                    &vxcfuncs::PkzbFuncGradientA,
                                    &vxcfuncs::PkzbFuncGradientB,
                                    &vxcfuncs::PkzbFuncHessianAB,
                                    &vxcfuncs::PkzbFuncHessianA,
                                    &vxcfuncs::PkzbFuncHessianB,
                                    &vxcfuncs::PkzbFuncCubicHessianAB,
                                    &vxcfuncs::PkzbFuncCubicHessianA,
                                    &vxcfuncs::PkzbFuncCubicHessianB);
    }

    void
    PkzbFuncGradientAB(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid)
    {
        // functional prefactors
        // determine number of grid points

        auto ngpoints = densityGrid.getNumberOfGridPoints();

        // set up pointers to density grid data

        auto rhoa = densityGrid.alphaDensity(0);

        auto rhob = densityGrid.betaDensity(0);

        auto grada = densityGrid.alphaDensityGradient(0);

        auto gradb = densityGrid.betaDensityGradient(0);

        auto taua = densityGrid.alphaDensitytau(0);

        auto taub = densityGrid.betaDensitytau(0);

        // set up pointers to functional data

        auto fexc = xcGradientGrid.xcFunctionalValues();

        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);

        auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);

        auto ggradb = xcGradientGrid.xcGradientValues(xcvars::gradb);

        auto gtaua = xcGradientGrid.xcGradientValues(xcvars::taua);

        auto gtaub = xcGradientGrid.xcGradientValues(xcvars::taub);

        #pragma omp simd aligned(rhoa, rhob, grada, gradb, taua, taub, fexc, grhoa, grhob, ggrada, ggradb, gtaua, gtaub : VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            // Constants

            const double kappa = 0.804;

            const double D = 0.113;

            const double c1 = 10.0 / 81.0;

            const double c2 = 146.0 / 2025.0;

            const double c3 = -73.0 / 405.0;

            const double c4 = D + 1.0 / kappa * std::pow(10.0 / 81.0, 2.0);

            double grada2 = grada[i] * grada[i];

            double gradb2 = gradb[i] * gradb[i];

            double pi = mathconst::getPiValue();

            double pi_sqr = std::pow(pi, 2.0);

            double C = -3.0 / (4.0 * pi) * std::pow(3.0 * pi_sqr, 1.0 / 3.0) * std::pow(2.0, 4.0 / 3.0);

            double eps_slater_a = C * std::pow(rhoa[i], 4.0 / 3.0);

            double eps_slater_b = C * std::pow(rhob[i], 4.0 / 3.0);

            double denom1 = 4.0 * std::pow(3.0 * pi_sqr, 2.0 / 3.0) * std::pow(2.0, 8.0 / 3.0);

            double Gamma_a = 4.0 * grada2 / (denom1 * std::pow(rhoa[i], 8.0 / 3.0));

            double Gamma_b = 4.0 * gradb2 / (denom1 * std::pow(rhob[i], 8.0 / 3.0));

            double denom2 = 2.0 * std::pow(3.0 * pi_sqr, 2.0 / 3.0) * std::pow(2.0, 5.0 / 3.0);

            double Upsilon_a = 3.0 * 2.0 * taua[i] / (denom2 * std::pow(rhoa[i], 5.0 / 3.0)) - 9.0 / 20.0 - Gamma_a / 12.0;

            double Upsilon_b = 3.0 * 2.0 * taub[i] / (denom2 * std::pow(rhob[i], 5.0 / 3.0)) - 9.0 / 20.0 - Gamma_b / 12.0;

            double Gamma_a2 = Gamma_a * Gamma_a;

            double Gamma_b2 = Gamma_b * Gamma_b;

            double Upsilon_a2 = Upsilon_a * Upsilon_a;

            double Upsilon_b2 = Upsilon_b * Upsilon_b;

            double Omega_a = c1 * Gamma_a + c2 * Upsilon_a2 + c3 * Upsilon_a * Gamma_a + c4 * Gamma_a2;

            double Omega_b = c1 * Gamma_b + c2 * Upsilon_b2 + c3 * Upsilon_b * Gamma_b + c4 * Gamma_b2;

            double F_a = 1 + kappa - kappa / (1.0 + Omega_a / kappa);

            double F_b = 1 + kappa - kappa / (1.0 + Omega_b / kappa);

            double eps_Pkzb_a = F_a * eps_slater_a;

            double eps_Pkzb_b = F_b * eps_slater_b;

            double C2 = -std::pow(3.0 * pi_sqr, 1 / 3.0) / pi * std::pow(2.0, 4.0 / 3.0);

            double diff_eps_a = C2 * std::pow(rhoa[i], 1.0 / 3.0);

            double diff_eps_b = C2 * std::pow(rhob[i], 1.0 / 3.0);

            double denom3 = 4.0 * std::pow(3.0 * pi_sqr, 2.0 / 3.0) * std::pow(2.0, 8.0 / 3.0);

            double diff_p_rho_a = -8.0 / 3.0 * 4.0 * grada2 / (denom3 * std::pow(rhoa[i], 11.0 / 3.0));

            double diff_p_rho_b = -8.0 / 3.0 * 4.0 * gradb2 / (denom3 * std::pow(rhob[i], 11.0 / 3.0));

            double denom4 = 2.0 * std::pow(3.0 * pi_sqr, 2.0 / 3.0) * std::pow(2.0, 5.0 / 3.0);

            double diff_q_rho_a = -5.0 / 3.0 * 6.0 * taua[i] / (denom4 * std::pow(rhoa[i], 8.0 / 3.0)) - 1 / 12.0 * diff_p_rho_a;

            double diff_q_rho_b = -5.0 / 3.0 * 6.0 * taub[i] / (denom4 * std::pow(rhob[i], 8.0 / 3.0)) - 1 / 12.0 * diff_p_rho_b;

            double diff_p_nabla_a = 8.0 * grada[i] / (denom3 * std::pow(rhoa[i], 8.0 / 3.0));

            double diff_p_nabla_b = 8.0 * gradb[i] / (denom3 * std::pow(rhob[i], 8.0 / 3.0));

            double diff_q_nabla_a = -1 / 12.0 * diff_p_nabla_a;

            double diff_q_nabla_b = -1 / 12.0 * diff_p_nabla_b;

            double diff_q_tau_a = 3.0 * 2.0 / (denom4 * std::pow(rhoa[i], 5.0 / 3.0));

            double diff_q_tau_b = 3.0 * 2.0 / (denom4 * std::pow(rhob[i], 5.0 / 3.0));

            double diff_x_rho_a = c1 * diff_p_rho_a + 2.0 * c2 * Upsilon_a * diff_q_rho_a + c3 * (diff_q_rho_a * Gamma_a + Upsilon_a * diff_p_rho_a) +
                                  2.0 * c4 * Gamma_a * diff_p_rho_a;

            double diff_x_rho_b = c1 * diff_p_rho_b + 2.0 * c2 * Upsilon_b * diff_q_rho_b + c3 * (diff_q_rho_a * Gamma_b + Upsilon_b * diff_p_rho_b) +
                                  2.0 * c4 * Gamma_b * diff_p_rho_b;

            double diff_x_nabla_a = c1 * diff_p_nabla_a + 2.0 * c2 * Upsilon_a * diff_q_nabla_a +
                                    c3 * (diff_q_nabla_a * Gamma_a + Upsilon_a * diff_p_nabla_a) + 2.0 * c4 * Gamma_a * diff_p_nabla_a;

            double diff_x_nabla_b = c1 * diff_p_nabla_b + 2.0 * c2 * Upsilon_b * diff_q_nabla_b +
                                    c3 * (diff_q_nabla_b * Gamma_b + Upsilon_b * diff_p_nabla_b) + 2.0 * c4 * Gamma_b * diff_p_nabla_b;

            double diff_x_tau_a = 2.0 * c2 * Upsilon_a * diff_q_tau_a + c3 * diff_q_tau_a * Gamma_a;

            double diff_x_tau_b = 2.0 * c2 * Upsilon_b * diff_q_tau_b + c3 * diff_q_tau_a * Gamma_b;

            double diff_F_x_a = std::pow(1 + Omega_a / kappa, -2);

            double diff_F_x_b = std::pow(1 + Omega_b / kappa, -2);

            double diff_F_rho_a = diff_F_x_a * diff_x_rho_a;

            double diff_F_rho_b = diff_F_x_b * diff_x_rho_b;

            double diff_F_nabla_a = diff_F_x_a * diff_x_nabla_a;

            double diff_F_nabla_b = diff_F_x_b * diff_x_nabla_b;

            double diff_F_tau_a = diff_F_x_a * diff_x_tau_a;

            double diff_F_tau_b = diff_F_x_b * diff_x_tau_b;

            fexc[i] += 0.5 * (eps_Pkzb_a + eps_Pkzb_b);

            grhoa[i] += 0.5 * (diff_eps_a * F_a + eps_slater_a * diff_F_rho_a);

            grhob[i] += 0.5 * (diff_eps_b * F_b + eps_slater_b * diff_F_rho_b);

            ggrada[i] += 0.5 * eps_slater_a * diff_F_nabla_a;

            ggradb[i] += 0.5 * eps_slater_b * diff_F_nabla_b;

            gtaua[i] += 0.5 * eps_slater_a * diff_F_tau_a;

            gtaub[i] += 0.5 * eps_slater_b * diff_F_tau_b;
        }
    }

    void
    PkzbFuncGradientA(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncGradientB(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncHessianAB(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncHessianA(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncHessianB(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncCubicHessianAB(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncCubicHessianA(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

    void
    PkzbFuncCubicHessianB(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid)
    {
    }

}  // namespace vxcfuncs
