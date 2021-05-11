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

#include "LYPFunctional.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    CXCFunctional
    setLYPFunctional()
    {
        return CXCFunctional({"LYP"}, xcfun::gga, 0.0, {setPrimitiveLYPFunctional()}, {1.0});
    }
    
    CPrimitiveFunctional
    setPrimitiveLYPFunctional()
    {
        return CPrimitiveFunctional({"LYP"}, xcfun::gga,
                                    &vxcfuncs::LYPFuncGradientAB,
                                    &vxcfuncs::LYPFuncGradientA,
                                    &vxcfuncs::LYPFuncGradientB,
                                    &vxcfuncs::LYPFuncHessianAB,
                                    &vxcfuncs::LYPFuncHessianA,
                                    &vxcfuncs::LYPFuncHessianB);
    }
    
    void
    LYPFuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                      const double           factor,
                      const CDensityGrid&    densityGrid)
    {
        double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
        
        double cf = 0.3 * std::pow(3 * mathconst::getPiValue() * mathconst::getPiValue(), 2.0 / 3.0);
        
        // determine number of grid points

        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        auto gradab = densityGrid.mixedDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);
        
        auto ggradb = xcGradientGrid.xcGradientValues(xcvars::gradb);
        
        auto ggradab = xcGradientGrid.xcGradientValues(xcvars::gradab);
        
        #pragma omp simd aligned(rhoa, rhob, grada, gradb, gradab, fexc, grhoa, grhob, ggrada, ggradb, ggradab: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i] + rhob[i];
            
            double rho2 = rho * rho;
            
            double rhoa2 = rhoa[i] * rhoa[i];
            
            double rhob2 = rhob[i] * rhob[i];
        
            double ngrada2 = grada[i] * grada[i];
            
            double ngradb2 = gradb[i] * gradb[i];
            
            double ngrad2 = ngrada2 + ngradb2 + 2.0 * gradab[i];
            
            double rhom13 = std::pow(rho, -1.0 / 3.0);
            
            double rho13 = std::pow(rho, 1.0 / 3.0);
            
            double drho13 = d + rho13;
            
            double drho3_2 = drho13 * drho13;
            
            double expcr = std::exp(-c / rho13);
            
            double grad = grada[i] + gradb[i];
            
            double grada2 = ngrada2;
            
            double gradb2 = ngradb2;
            
            double grad2 = ngrad2;
            
            double sA = (rhoa[i] * grada2 + rhob[i] * gradb2) / rho;
            
            double sA10 = (grada[i] - gradb[i]) * grad * rhob[i] / rho2;
            
            double sA01 = (gradb[i] - grada[i]) * grad * rhoa[i] / rho2;
            
            double denom = 1.0 + d * rhom13;

            double omega = std::exp(-c * rhom13) / denom * pow(rho, -11.0 / 3.0);
    
            double delta = rhom13 * (c + d / denom);
        
            double t1 = std::pow(2.0, 11.0 / 3.0) * cf * (std::pow(rhoa[i], 8.0 / 3.0) + std::pow(rhob[i], 8.0 / 3.0));
            
            double t2 = (47.0 - 7.0 * delta) * ngrad2 / 18.0;
        
            double t3 = -(2.5 - delta / 18.0) * (ngrada2 + ngradb2);
            
            double t4 = (11.0 - delta) / 9.0 * (rhoa[i] * ngrada2 + rhob[i] * ngradb2) / rho;
            
            double t5 = -2.0 / 3.0 * rho2 * ngrad2;
            
            double t6 = ((2.0 / 3.0 * rho2 - rhoa[i] * rhoa[i]) * ngradb2 + (2.0 / 3.0 * rho2 - rhob[i] * rhob[i]) * ngrada2);
            
            fexc[i] += -factor * a * (4.0 * rhoa[i] * rhob[i] / (denom * rho) + b * omega * (rhoa[i] * rhob[i] * (t1 + t2 + t3 + t4) + t5 + t6));
            
            double om = expcr * std::pow(rho, -11.0 / 3.0) / (1.0 +d / rho13);
            
            double om_1 = (-11.0 * rho * rho13 + rho * (c - 10.0 * d) + rho / rho13 * c * d) / (3.0 * std::pow(rho, 16.0 / 3.0) * drho3_2) * expcr;
            
            double dl = c / rho13 + d / (rho13 * (1 + d / rho13));
            
            double dl_1 = (-c * drho3_2 - d * rho / rho13) / (3.0 * rho * rho13 * drho3_2);
            
            double f0_1000 = 4.0 * rhob[i] * (d * rhoa[i] + 3.0 * rhob[i] * drho13) * rho13 / (3.0 * rho * rho * drho3_2);
            
            double f0_0100 = 4.0 * rhoa[i] * (d * rhob[i] + 3.0 * rhoa[i] * drho13) * rho13 / (3.0 * rho * rho * drho3_2);
            
            double f1 = std::pow(2.0, 11.0 / 3.0) * cf * (std::pow(rhoa[i], 8.0 / 3.0) + std::pow(rhob[i], 8.0 / 3.0))
            
                      + (47.0 - 7.0 * dl) *grad2 / 18.0 + (dl - 45.0) * (grada2 + gradb2) / 18.0 + (11.0 - dl) * sA / 9.0;
            
            double f1_1000 = std::pow(2.0, 11.0 / 3.0) * cf * 8.0 / 3.0 * std::pow(rhoa[i], 5.0 / 3.0)
            
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_1 / 18.0 + (11.0 - dl) * sA10 / 9.0;
            
            double f1_0100 = std::pow(2.0, 11.0 / 3.0) * cf * 8.0 / 3.0 * std::pow(rhob[i], 5.0 / 3.0)
            
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_1 / 18.0 + (11.0 - dl) * sA01 / 9.0;
            
            double f1_0010 = (47.0 - 7.0 * dl) * grada[i] / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhoa[i] / rho) * grada[i] / 9.0;
            
            double f1_0001 = (47.0 - 7.0 * dl) * gradb[i] / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhob[i] / rho) * gradb[i] / 9.0;
            
            double f2 = -2.0 / 3.0 * rho2 * grad2 + (2.0 / 3.0 * rho2 - rhoa2) * gradb2 + (2.0 / 3.0 * rho2 - rhob2) * grada2;
            
            double f2_1000 = -8.0 / 3.0 * rho * gradab[i] - 2.0 * rhoa[i] * gradb2;
            
            double f2_0100 = -8.0 / 3.0 * rho * gradab[i] - 2.0 * rhob[i] * grada2;
            
            double f2_0010 = -2.0 * rhob2 * grada[i];
            
            double f2_0001 = -2.0 * rhoa2 * gradb[i];
            
            double f1_00001 = (47.0 - 7.0 * dl) / 9.0;
            
            double f2_00001 = -4.0 / 3.0 * rho2;
            
            grhoa[i] += factor * (-a * f0_1000 - a * b * om_1 * (rhoa[i] * rhob[i] * f1 + f2) - a * b * om * (rhob[i] * f1 + rhoa[i] * rhob[i] * f1_1000 + f2_1000));
            
            grhob[i] += factor * (-a * f0_0100 - a * b * om_1 * (rhoa[i] * rhob[i] * f1 + f2) - a * b * om * (rhoa[i] * f1 + rhoa[i] * rhob[i] * f1_0100 + f2_0100));
            
            ggrada[i] += factor * (-a * b * om * (rhoa[i] * rhob[i] * f1_0010 + f2_0010));
            
            ggradb[i] += factor * (-a * b * om * (rhoa[i] * rhob[i] * f1_0001 + f2_0001));
            
            ggradab[i] += factor * (-a * b * om * (rhoa[i] * rhob[i] * f1_00001 + f2_00001));
        }
    }
    
    void
    LYPFuncGradientA(      CXCGradientGrid& xcGradientGrid,
                     const double           factor,
                     const CDensityGrid&    densityGrid)
    {

        double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
        
        double cf = 0.3 * std::pow(3 * mathconst::getPiValue() * mathconst::getPiValue(), 2.0 / 3.0);
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        //auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        //auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        //auto gradab = densityGrid.mixedDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
                        
        auto ggradab = xcGradientGrid.xcGradientValues(xcvars::gradab);
        
        #pragma omp simd aligned(rhob, gradb, fexc, grhoa, grhob, ggradab: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhob[i];
            
            double rho2 = rho * rho;
                        
            double ngradb2 = gradb[i] * gradb[i];
            
            double ngrad2 =  ngradb2 ;
            
            double rhom13 = std::pow(rho, -1.0 / 3.0);
            
            double rho13 = std::pow(rho, 1.0 / 3.0);
            
            double drho13 = d + rho13;
            
            double drho3_2 = drho13 * drho13;
            
            double expcr = std::exp(-c / rho13);
                        
            double gradb2 = ngradb2;
            
            double grad2 = ngrad2;
            
            double sA = (rhob[i] * gradb2) / rho;
                        
            double denom = 1.0 + d * rhom13;

            double omega = std::exp(-c * rhom13) / denom * pow(rho, -11.0 / 3.0);
            
            double t5 = -2.0 / 3.0 * rho2 * ngrad2;
            
            double t6 = ((2.0 / 3.0 * rho2 ) * ngradb2 );
            
            fexc[i] += -factor * a * (b * omega *  (t5 + t6));
            
            double om = expcr * std::pow(rho, -11.0 / 3.0) / (1.0 +d / rho13);
            
            double om_1 = (-11.0 * rho * rho13 + rho * (c - 10.0 * d) + rho / rho13 * c * d) / (3.0 * std::pow(rho, 16.0 / 3.0) * drho3_2) * expcr;
            
            double dl = c / rho13 + d / (rho13 * (1 + d / rho13));
                        
            double f0_1000 = 4.0 * rhob[i] * (3.0 * rhob[i] * drho13) * rho13 / (3.0 * rho * rho * drho3_2);
            
            double f1 = std::pow(2.0, 11.0 / 3.0) * cf * (std::pow(rhob[i], 8.0 / 3.0))
            
                      + (47.0 - 7.0 * dl) *grad2 / 18.0 + (dl - 45.0) * (gradb2) / 18.0 + (11.0 - dl) * sA / 9.0;

            double f2 = -2.0 / 3.0 * rho2 * grad2 + (2.0 / 3.0 * rho2 ) * gradb2 ;
                      
            double f2_00001 = -4.0 / 3.0 * rho2;
            
            grhoa[i] += factor * (-a * f0_1000 - a * b * om_1 * (f2) - a * b * om * (rhob[i] * f1 ));
            
            grhob[i] += factor * (- a * b * om_1 * f2);
                                
            ggradab[i] += factor * (-a * b * om * ( f2_00001));
        }
    }

    void
    LYPFuncGradientB(      CXCGradientGrid& xcGradientGrid,
                     const double           factor,
                     const CDensityGrid&    densityGrid)
    {
        double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
        
        double cf = 0.3 * std::pow(3 * mathconst::getPiValue() * mathconst::getPiValue(), 2.0 / 3.0);
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
                
        auto grada = densityGrid.alphaDensityGradient(0);
        
        // set up pointers to functional data

        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        auto ggradab = xcGradientGrid.xcGradientValues(xcvars::gradab);

        #pragma omp simd aligned(rhoa, grada, fexc, grhoa, grhob,ggradab: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i];
            
            double rho2 = rho * rho;
                                
            double ngrada2 = grada[i] * grada[i];
                        
            double ngrad2 = ngrada2;
            
            double rhom13 = std::pow(rho, -1.0 / 3.0);
            
            double rho13 = std::pow(rho, 1.0 / 3.0);
            
            double drho13 = d + rho13;
            
            double drho3_2 = drho13 * drho13;
            
            double expcr = std::exp(-c / rho13);
                        
            double grada2 = ngrada2;
                        
            double grad2 = ngrad2;
            
            double sA = (rhoa[i] * grada2) / rho;
                        
            double denom = 1.0 + d * rhom13;

            double omega = std::exp(-c * rhom13) / denom * pow(rho, -11.0 / 3.0);
                
            double t5 = -2.0 / 3.0 * rho2 * ngrad2;
            
            double t6 = ((2.0 / 3.0 * rho2) * ngrada2);
            
            fexc[i] += -factor * a * (b * omega * (t5 + t6));
            
            double om = expcr * std::pow(rho, -11.0 / 3.0) / (1.0 +d / rho13);
            
            double om_1 = (-11.0 * rho * rho13 + rho * (c - 10.0 * d) + rho / rho13 * c * d) / (3.0 * std::pow(rho, 16.0 / 3.0) * drho3_2) * expcr;
            
            double dl = c / rho13 + d / (rho13 * (1 + d / rho13));
                        
            double f0_0100 = 4.0 * rhoa[i] * (3.0 * rhoa[i] * drho13) * rho13 / (3.0 * rho * rho * drho3_2);
            
            double f1 = std::pow(2.0, 11.0 / 3.0) * cf * (std::pow(rhoa[i], 8.0 / 3.0)) + (47.0 - 7.0 * dl) *grad2 / 18.0 + (dl - 45.0) * (grada2) / 18.0 + (11.0 - dl) * sA / 9.0;
                        
            double f2 = -2.0 / 3.0 * rho2 * grad2 + (2.0 / 3.0 * rho2 ) * grada2;                                    
            
            double f2_00001 = -4.0 / 3.0 * rho2;
            
            grhoa[i] += factor * (- a * b * om_1 * (f2));
            
            grhob[i] += factor * (-a * f0_0100 - a * b * om_1 * (f2) - a * b * om * (rhoa[i] * f1 ));
                                    
            ggradab[i] += factor * (-a * b * om * (f2_00001));
        }
    }
    
    void
    LYPFuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                     const double          factor,
                     const CDensityGrid&   densityGrid)
    {
        double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
        
        double CF = 0.3 * std::pow(3.0 * mathconst::getPiValue() * mathconst::getPiValue(), 2.0 / 3.0);
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        auto gradab = densityGrid.mixedDensityGradient(0);
        
        // set up pointers to hessian data
        
        auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);
        
        auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        auto ggrad_aa = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);
        
        auto ggrad_bb = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradb);
        
        auto gmix_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);
        
        auto gmix_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradb);
        
        auto gmix_ba = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::grada);
        
        auto gmix_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradb);
        
        auto gmix_ac = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradab);
        
        auto gmix_bc = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradab);
        
        #pragma omp simd aligned(rhoa, rhob, grada, gradb, grho_aa, grho_ab, grho_bb, ggrad_aa, ggrad_bb,\
                                gmix_aa, gmix_ab, gmix_ba, gmix_bb, gmix_ac, gmix_bc: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i] + rhob[i];
            
            double rho2 = rho * rho;
            
            double rhoa2 = rhoa[i] * rhoa[i];
            
            double rhob2 = rhob[i] * rhob[i];
            
            double rhoab = rhoa[i] * rhob[i];
            
            double grad = grada[i] + gradb[i];
            
            double rho13  = std::pow(rho, 1.0 / 3.0);
            
            double drho13 = D + rho13;
            
            double drho3_2 = drho13 * drho13;
            
            double drho3_3 = drho13 * drho3_2;
            
            double expcr  = std::exp(-C / rho13);
            
            double de = C / ( 3.0 * rho * rho13);
            
            double de1 = -4.0 * C / ( 9.0 * rho2 * rho13);
            
            double grada2 = grada[i] * grada[i];
            
            double gradb2 = gradb[i] * gradb[i];
            
            double grad2 = grada2 + gradb2 + 2.0 * gradab[i];
           
            /* s derivatives */
            
            double sA = (rhoa[i] * grada2 + rhob[i] * gradb2) / rho;
            
            double sA1000 = (grada[i] - gradb[i]) * grad * rhob[i] / rho2;
            
            double sA0100 = (gradb[i] - grada[i]) * grad * rhoa[i] / rho2;
            
            double sA2000 =  -2.0 * sA1000 / rho;
            
            double sA0200 =  -2.0 * sA0100 / rho;
            
            double sA1100 = (rhoa[i] - rhob[i]) * (grada[i] - gradb[i]) * grad / (rho2 * rho);
            
            /* omega derivatives */
            
            double omx  = std::pow(rho, -11.0 / 3.0) / (1.0 + D / rho13);
            
            double omx1 = (-11.0 * rho13 - 10.0 * D) / (3.0 * rho2 * rho2 * rho13 * drho3_2);
            
            double omx2 = 2.0 * (77.0 * rho / rho13 + 141.0 * D * rho13 + 65.0 * D * D) / (9.0 * std::pow(rho, 16.0 / 3.0) * drho3_3);
            
            double om = expcr * omx;
            
            double om_1 = expcr * (omx1 + omx * de);
            
            double om_2 = expcr * (omx2 + 2.0 * omx1 * de + omx * (de1 + de * de));
            
            /* delta derivatives */
            
            double dl = C / rho13 + D / (rho13 * (1.0 + D / rho13));
            
            double dl_1 = (-C * drho3_2 - D * rho / rho13) / (3.0 * rho * rho13 * drho3_2);
            
            double dl_2 = (4.0 * rho * (C + D) + 2.0 * D * (6.0 * rho13 * C * D + 2.0 * C * D * D + rho / rho13 * (6 * C + D))) / (9.0 * rho2 * rho13 * drho3_3);
            
            /* f0 derivatives */
            
            double f0_2000 = -8.0 * rhob[i] * (rhoa[i] * D * (D + 2.0 * rho13) + 3.0 * rhob[i] * drho13 * (2.0 * D + 3.0 * rho13)) / (9.0 * rho2 * rho / rho13 * drho3_2 * drho13);
            
            double f0_0200 = -8.0 * rhoa[i] * (rhob[i] * D * (D + 2.0 * rho13) + 3.0 * rhoa[i] * drho13 * (2.0 * D + 3.0 * rho13)) / (9.0 * rho2 * rho / rho13 * drho3_2 * drho13);
            
            double f0_1100 = 4.0 * (18.0 * rhoa[i] * rhob[i] * rho / rho13 + rho13 * (3.0 * rhoa2 + 32.0 * rhoab + 3.0 * rhob2) * D
                                    
                           + (3.0 * rhoa2 + 16.0 * rhoab + 3.0 * rhob2) * D * D) / (9.0 * rho2 * rho / rho13 * drho3_3);
            
            /* f1 derivatives */
            
            double f1 = std::pow(2.0, 11.0 / 3.0) * CF * (std::pow(rhoa[i], 8.0 / 3.0) + std::pow(rhob[i], 8.0 / 3.0))
            
                      + (47.0 - 7.0 * dl) * grad2 / 18.0 + (dl - 45.0) * (grada2 + gradb2) / 18.0 + (11.0 - dl) * sA / 9.0;
            
            double f1_1000 = std::pow(2.0, 11.0 / 3.0) * CF * 8.0 / 3.0 * std::pow(rhoa[i], 5.0 / 3.0)
           
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_1 / 18.0 + (11.0 - dl) * sA1000 / 9.0;
            
            double f1_0100 = std::pow(2.0, 11.0 / 3.0) * CF * 8.0 / 3.0 * std::pow(rhob[i], 5.0 / 3.0)
           
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_1 / 18.0 + (11.0 - dl) * sA0100 / 9.0;
            
            double f1_0010 = (47.0 - 7.0 * dl) * grada[i] / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhoa[i] / rho) * grada[i] / 9.0;
            
            double f1_0001 = (47.0 - 7.0 * dl) * gradb[i] / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhob[i] / rho) * gradb[i] / 9.0;
            
            double f1_00001 = (47.0 - 7.0 * dl) / 9.0;
            
            double f1_2000 = std::pow(2.0, 11.0 / 3.0) * CF * 40.0 / 9.0 * std::pow(rhoa[i], 2.0 / 3.0) - 2.0 * sA1000 * dl_1 / 9.0
            
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_2 / 18.0 + (11.0 - dl) * sA2000 / 9.0;
            
            double f1_0200 = std::pow(2.0, 11.0 / 3.0) * CF * 40.0 / 9.0 * std::pow(rhob[i], 2.0 / 3.0) - 2.0 * sA0100 * dl_1 / 9.0
            
                           + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_2 / 18.0 + (11.0 - dl) * sA0200 / 9.0;
            
            double f1_0020 = (47.0 - 7.0 * dl) / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhoa[i] / rho) / 9.0;
            
            double f1_0002 = (47.0 - 7.0 * dl) / 9.0 + (dl - 45.0 + (22.0 - 2.0 * dl) * rhob[i] / rho) / 9.0;
            
            double f1_1100 = -2.0 * sA0100 * dl_1 / 18.0 + (grada2 + gradb2 - 7.0 * grad2 - 2.0 * sA) * dl_2 / 18.0
            
                           - dl_1 * sA1000 / 9.0 + (11.0 - dl) * sA1100 / 9.0;
            
            double f1_1010 = (grada[i] * (1.0 - 2.0 * rhoa[i] / rho) - 7.0 * grada[i]) * dl_1 / 9.0 + (11.0 - dl) * rhob[i] / rho2 * grada[i] / 4.5;
            
            double f1_0101 = (gradb[i] * (1.0 - 2.0 * rhob[i] / rho) - 7.0 * gradb[i]) * dl_1 / 9.0 + (11.0 - dl) * rhoa[i] / rho2 * gradb[i] / 4.5;
            
            double f1_1001 = (gradb[i] * (1.0 - 2.0 * rhob[i] / rho) - 7.0 * gradb[i]) * dl_1 / 9.0 - (11.0 - dl) * rhob[i] / rho2 * gradb[i] / 4.5;
            
            double f1_0110 = (grada[i] * (1.0 - 2.0 * rhoa[i] / rho) - 7.0 * grada[i]) * dl_1 / 9.0 - (11.0 - dl) * rhoa[i] / rho2 * grada[i] / 4.5;
            
            double f1_10001 = -7.0 * dl_1 / 9.0;
            
            double f1_01001 = -7.0 * dl_1 / 9.0;
            
            /* f2 derivatives */
            
            double f2 = -2.0 / 3.0 * rho2  * grad2 + (2.0 / 3.0 * rho2 - rhoa2) * gradb2 + (2.0 / 3.0 * rho2 - rhob2) * grada2;
            
            double f2_1000 = -8.0 / 3.0 * rho * gradab[i] - 2.0 * rhoa[i] * gradb2;
            
            double f2_0100 = -8.0 / 3.0 * rho * gradab[i] - 2.0 * rhob[i] * grada2;
            
            double f2_0010 = -2.0 * rhob2 * grada[i];
            
            double f2_0001 = -2.0 * rhoa2 * gradb[i];
            
            double f2_00001 = -4.0 / 3.0 * rho2;
            
            double f2_2000 = -8.0 / 3.0 * gradab[i] - 2.0 * gradb2;
            
            double f2_0200 = -8.0 / 3.0 * gradab[i] - 2.0 * grada2;
            
            double f2_0020 = -2.0 * rhob2;
            
            double f2_0002 = -2.0 * rhoa2;
            
            double f2_1100  = -8.0 / 3.0 * gradab[i];
            
            double f2_1010 = 0.0;
            
            double f2_0101 = 0.0;
            
            double f2_1001 = -4.0 * rhoa[i] * gradb[i];
            
            double f2_0110 = -4.0 * rhob[i] * grada[i];
            
            double f2_10001 = -8.0 / 3.0 * rho;
            
            double f2_01001 = -8.0 / 3.0 * rho;
            
            /* derivatives sums */
            
            double rff = rhoa[i] * rhob[i] * f1 + f2;
            
            double rff_1000 = rhob[i] * f1 + rhoa[i] * rhob[i] * f1_1000 + f2_1000;
            
            double rff_0100 = rhoa[i] * f1 + rhoa[i] * rhob[i] * f1_0100 + f2_0100;
            
            double rff_2000 = 2.0 * rhob[i] * f1_1000 + rhoa[i] * rhob[i] * f1_2000 + f2_2000;
            
            double rff_0200 = 2.0 * rhoa[i] * f1_0100 + rhoa[i] * rhob[i] * f1_0200 + f2_0200;
            
            double rff_1100 = f1 + rhob[i] * f1_0100 + rhoa[i] * f1_1000 + rhoa[i] * rhob[i] * f1_1100 + f2_1100;
            
            double rff_0010 = rhoa[i] * rhob[i] * f1_0010 + f2_0010;
            
            double rff_0001 = rhoa[i] * rhob[i] * f1_0001 + f2_0001;
            
            double rff_1010 = rhob[i] * f1_0010 + rhoa[i] * rhob[i] * f1_1010 + f2_1010;
            
            double rff_0101 = rhoa[i] * f1_0001 + rhoa[i] * rhob[i] * f1_0101 + f2_0101;
            
            /* derivatives sum with respect grada*gradb */
            
            double rff_00001 = rhoa[i] * rhob[i] * f1_00001 + f2_00001;
            
            double rff_10001 = rhob[i] * f1_00001 + rhoa[i] * rhob[i] * f1_10001 + f2_10001;
            
            double rff_01001 = rhoa[i] * f1_00001 + rhoa[i] * rhob[i] * f1_01001 + f2_01001;
            
            /* the final section: second derivatives */
            
            grho_aa[i] += factor * (-A * f0_2000 - A * B * (om_2 * rff + 2.0 * om_1 * rff_1000 + om * rff_2000));

            grho_bb[i] += factor * (-A * f0_0200 - A * B * (om_2 * rff + 2.0 * om_1 * rff_0100 + om * rff_0200));
            
            ggrad_aa[i] += factor * (-A * B * om * (rhoa[i] * rhob[i] * f1_0020 + f2_0020));
            
            ggrad_bb[i] += factor * (-A * B * om * (rhoa[i] * rhob[i] * f1_0002 + f2_0002));
            
            /* the mixed derivatives */
            
            grho_ab[i] += factor*(-A * f0_1100 - A * B * (om_2 * rff + om_1 * rff_0100 + om_1 * rff_1000 + om * rff_1100));
            
            gmix_aa[i] += factor * (-A * B * (om_1 * rff_0010 + om * rff_1010));
            
            gmix_ab[i] += factor * (-A * B * (om_1 * (rhoa[i] * rhob[i] * f1_0001 + f2_0001) + om * (rhob[i] * f1_0001 + rhoa[i] * rhob[i] * f1_1001 + f2_1001)));
            
            gmix_bb[i] += factor * (-A * B * (om_1 * rff_0001 + om * rff_0101));
            
            gmix_ba[i] += factor * (-A * B * (om_1 * (rhoa[i] * rhob[i] * f1_0010 + f2_0010) + om * (rhoa[i] * f1_0010 + rhoa[i] * rhob[i] * f1_0110 + f2_0110)));
            
            gmix_ac[i] += factor * (-A * B * (om_1 * rff_00001 + om * rff_10001));
            
            gmix_bc[i] += factor * (-A * B * (om_1 * rff_00001 + om * rff_01001));
        }
    }
    
    void
    LYPFuncHessianA(      CXCHessianGrid& xcHessianGrid,
                    const double          factor,
                    const CDensityGrid&   densityGrid)
    {
        
    }
    
    void
    LYPFuncHessianB(      CXCHessianGrid& xcHessianGrid,
                    const double          factor,
                    const CDensityGrid&   densityGrid)
    {
        
    }
    
}  // namespace vxcfuncs
