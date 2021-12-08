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

#include "VWN3Functional.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    CXCFunctional
    setVWN3Functional()
    {
        return CXCFunctional({"VWN3"}, xcfun::lda, 0.0, {setPrimitiveVWN3Functional()}, {1.0});
    }
    
    CPrimitiveFunctional
    setPrimitiveVWN3Functional()
    {
        return CPrimitiveFunctional({"VWN3"}, xcfun::lda,
                                    &vxcfuncs::VWN3FuncGradientAB,
                                    &vxcfuncs::VWN3FuncGradientA,
                                    &vxcfuncs::VWN3FuncGradientB,
                                    &vxcfuncs::VWN3FuncHessianAB,
                                    &vxcfuncs::VWN3FuncHessianA,
                                    &vxcfuncs::VWN3FuncHessianB,
                                    &vxcfuncs::VWN3FuncCubicHessianAB,
                                    &vxcfuncs::VWN3FuncCubicHessianA,
                                    &vxcfuncs::VWN3FuncCubicHessianB);
    }
    
    void
    VWN3FuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                       const double           factor,
                       const CDensityGrid&    densityGrid)
    {
        
        const double spinpolf = 1.92366105093154;

        const double fourthree   = 1.333333333333333;
        
        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);

        // Ferromagnetic fitting parameters

        double pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578, px0_f = -0.7432940;
        
        double pq_f = std::sqrt(4.0 * pc_f - pb_f* pb_f);
        
        double pxf0_f = px0_f * px0_f + pb_f * px0_f + pc_f;
        
        double pyf0_f = pq_f / (pb_f + 2.0 * px0_f);
        
        double b_f = px0_f / pxf0_f;
        
        double c_f = pxf0_f * pyf0_f;
        
        double acon_f = b_f * pb_f - 1.0;
        
        double bcon_f = 2.0 * acon_f + 2.0;
        
        double ccon_f = 2.0 * pb_f * (1.0 / pq_f - px0_f / c_f);

        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);

        auto rhob = densityGrid.betaDensity(0);

        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);

        // diamagnetic contribution
        
        #pragma omp simd aligned(rhoa, rhob, fexc, grhoa, grhob: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i] + rhob[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);

            double zeta = (rhoa[i]-rhob[i])/rho;

            double f_zeta  = spinpolf * (std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf * 4.0/3.0 * (std::pow(1+zeta,1.0/3.0)-std::pow(1-zeta,1.0/3.0));

            // Paramagnetic

            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpe1 = 2.0 * std::log(x) + acon * std::log(xf) - bcon * std::log(x - px0) + ccon * std::atan(yf);

            double eps_c0 = fpe1*0.5*pa;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);

            double diff_eps_c0 = 0.5*pa*(fpe1 + rho*fpex1*xrho);

            // Ferromagnetic

            double xf_f = x * x + pb_f * x + pc_f;
            
            double xfx_f = 2.0 * x + pb_f;
            
            double yf_f = pq_f / xfx_f;
            
            double fpe1_f = 2.0 * std::log(x) + acon_f * std::log(xf_f) - bcon_f * std::log(x - px0_f) + ccon_f * std::atan(yf_f);

            double eps_c1 = fpe1_f*0.5*pa_f;
            
            double fpex1_f = 2.0 / x + acon_f * xfx_f / xf_f - bcon_f / (x - px0) - ccon_f * (2.0 * yf_f / xfx_f) / (1.0 + yf_f * yf_f);
            
            double diff_eps_c1 = 0.5*pa_f*(fpe1_f + rho*fpex1_f*xrho);

            double vcfp = f_zeta*(diff_eps_c1-diff_eps_c0);

            double delta = f_zet1*(eps_c1-eps_c0);

            fexc[i] += (eps_c0 + f_zeta*(eps_c1-eps_c0))*rho*factor ;
            
            grhoa[i] +=  (diff_eps_c0 + vcfp + delta*(1-zeta))*factor;

            grhob[i] +=  (diff_eps_c0 + vcfp - delta*(1+zeta))*factor;

        }
    }
    
    void
    VWN3FuncGradientA(      CXCGradientGrid& xcGradientGrid,
                      const double           factor,
                      const CDensityGrid&    densityGrid)
    {

        const double spinpolf = 1.92366105093154;

        const double fourthree   = 1.333333333333333;

        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);
        
        // Ferromagnetic fitting parameters

        double pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578, px0_f = -0.7432940;
        
        double pq_f = std::sqrt(4.0 * pc_f - pb_f* pb_f);
        
        double pxf0_f = px0_f * px0_f + pb_f * px0_f + pc_f;
        
        double pyf0_f = pq_f / (pb_f + 2.0 * px0_f);
        
        double b_f = px0_f / pxf0_f;
        
        double c_f = pxf0_f * pyf0_f;
        
        double acon_f= b_f * pb_f - 1.0;
        
        double bcon_f = 2.0 * acon_f + 2.0;
        
        double ccon_f = 2.0 * pb_f * (1.0 / pq_f - px0_f / c_f);

        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
                
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);

        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

        // diamagnetic contribution
        
        #pragma omp simd aligned(rhob, fexc, grhoa, grhob: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhob[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);

            double zeta = -1;

            double f_zeta  = spinpolf*(std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf*4.0/3.0*(std::pow(1+zeta,1.0/3.0)-std::pow(1-zeta,1.0/3.0));

            // Paramagnetic

            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpe1 = 2.0 * std::log(x) + acon * std::log(xf) - bcon * std::log(x - px0) + ccon * std::atan(yf);

            double eps_c0 = fpe1*0.5*pa;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);

            double diff_eps_c0 = 0.5*pa*(fpe1 + rho*fpex1*xrho);

            // Ferromagnetic

            double xf_f = x * x + pb_f * x + pc_f;
            
            double xfx_f = 2.0 * x + pb_f;
            
            double yf_f = pq_f / xfx_f;
            
            double fpe1_f = 2.0 * std::log(x) + acon_f * std::log(xf_f) - bcon_f * std::log(x - px0_f) + ccon_f * std::atan(yf_f);

            double eps_c1 = fpe1_f*0.5*pa_f;
            
            double fpex1_f = 2.0 / x + acon_f * xfx_f / xf_f - bcon_f / (x - px0) - ccon_f * (2.0 * yf_f / xfx_f) / (1.0 + yf_f * yf_f);
            
            double diff_eps_c1 = 0.5*pa_f*(fpe1_f + rho*fpex1_f*xrho);

            double vcfp = f_zeta*(diff_eps_c1-diff_eps_c0);

            double delta = f_zet1*(eps_c1-eps_c0);

            fexc[i] += (eps_c0 + f_zeta*(eps_c1-eps_c0))*rho*factor ;
            
            grhoa[i] +=  (diff_eps_c0 + vcfp + delta*(1-zeta))*factor;

            grhob[i] +=  (diff_eps_c0 + vcfp - delta*(1+zeta))*factor;
        }
    }
    
    void
    VWN3FuncGradientB(      CXCGradientGrid& xcGradientGrid,
                      const double           factor,
                      const CDensityGrid&    densityGrid)
    {


        const double spinpolf = 1.92366105093154;

        const double fourthree   = 1.333333333333333;

        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);
        
        // Ferromagnetic fitting parameters

        double pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578, px0_f = -0.7432940;
        
        double pq_f = std::sqrt(4.0 * pc_f - pb_f* pb_f);
        
        double pxf0_f = px0_f * px0_f + pb_f * px0_f + pc_f;
        
        double pyf0_f = pq_f / (pb_f + 2.0 * px0_f);
        
        double b_f = px0_f / pxf0_f;
        
        double c_f = pxf0_f * pyf0_f;
        
        double acon_f= b_f * pb_f - 1.0;
        
        double bcon_f = 2.0 * acon_f + 2.0;
        
        double ccon_f = 2.0 * pb_f * (1.0 / pq_f - px0_f / c_f);

        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
                
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);

        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);

        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        // diamagnetic contribution
        
        #pragma omp simd aligned(rhoa,fexc, grhoa, grhob: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);

            double zeta = 1;

            double f_zeta  = spinpolf*(std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf*4.0/3.0*(std::pow(1+zeta,1.0/3.0)-std::pow(1-zeta,1.0/3.0));

            // Paramagnetic

            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpe1 = 2.0 * std::log(x) + acon * std::log(xf) - bcon * std::log(x - px0) + ccon * std::atan(yf);

            double eps_c0 = fpe1*0.5*pa;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);

            double diff_eps_c0 = 0.5*pa*(fpe1 + rho*fpex1*xrho);

            // Ferromagnetic

            double xf_f = x * x + pb_f * x + pc_f;
            
            double xfx_f = 2.0 * x + pb_f;
            
            double yf_f = pq_f / xfx_f;
            
            double fpe1_f = 2.0 * std::log(x) + acon_f * std::log(xf_f) - bcon_f * std::log(x - px0_f) + ccon_f * std::atan(yf_f);

            double eps_c1 = fpe1_f*0.5*pa_f;
            
            double fpex1_f = 2.0 / x + acon_f * xfx_f / xf_f - bcon_f / (x - px0) - ccon_f * (2.0 * yf_f / xfx_f) / (1.0 + yf_f * yf_f);
            
            double diff_eps_c1 = 0.5*pa_f*(fpe1_f + rho*fpex1_f*xrho);

            double vcfp = f_zeta*(diff_eps_c1-diff_eps_c0);

            double delta = f_zet1*(eps_c1-eps_c0);

            fexc[i] += (eps_c0 + f_zeta*(eps_c1-eps_c0))*rho*factor ;
            
            grhoa[i] +=  (diff_eps_c0 + vcfp + delta*(1-zeta))*factor;

            grhob[i] +=  (diff_eps_c0 + vcfp - delta*(1+zeta))*factor;

        }
    }
    
    void
    VWN3FuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {
        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;

        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);
        
        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
        
        double fpre = factor * 0.5 * pa;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);
        
        auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        // diamagnetic contribution
        
        #pragma omp simd aligned(rhoa, rhob, grho_aa, grho_ab, grho_bb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i] + rhob[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);
            
            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);
            
            double fpexx1 = -2.0 / (x * x) + acon * (2.0 / xf - xfx * xfx / (xf * xf))
           
                          + bcon / ((x - px0) * (x - px0)) + ccon * 8.0 * pq * xfx / ((pq * pq + xfx * xfx) * (pq * pq + xfx * xfx));
            
            double gval = fpre * xrho * (2.0 * fpex1 + rho * fpexx1 * xrho - fpex1 * 7.0 / 6.0);
            
            grho_aa[i] += gval;
            
            grho_ab[i] += gval;
            
            grho_bb[i] += gval;
        }
    }
    
    void
    VWN3FuncHessianA(      CXCHessianGrid& xcHessianGrid,
                     const double          factor,
                     const CDensityGrid&   densityGrid)
    {
        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);
        
        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
        
        double fpre = factor * 0.5 * pa;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        // diamagnetic contribution
        
        #pragma omp simd aligned(rhob, grho_bb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhob[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);
            
            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);
            
            double fpexx1 = -2.0 / (x * x) + acon * (2.0 / xf - xfx * xfx / (xf * xf))
           
                          + bcon / ((x - px0) * (x - px0)) + ccon * 8.0 * pq * xfx / ((pq * pq + xfx * xfx) * (pq * pq + xfx * xfx));
            
            double gval = fpre * xrho * (2.0 * fpex1 + rho * fpexx1 * xrho - fpex1 * 7.0 / 6.0);
            
            grho_bb[i] += gval;
        }
    }
    
    void
    VWN3FuncHessianB(      CXCHessianGrid& xcHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {
        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);
        
        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
        
        double fpre = factor * 0.5 * pa;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        // diamagnetic contribution
        
        #pragma omp simd aligned(rhoa, grho_aa: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i];
            
            double x = dcrs * std::pow(rho, f16);
            
            double xrho  = dcrs * f16 * std::pow(rho, f76);
            
            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;
            
            double yf = pq / xfx;
            
            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);
            
            double fpexx1 = -2.0 / (x * x) + acon * (2.0 / xf - xfx * xfx / (xf * xf))
           
                          + bcon / ((x - px0) * (x - px0)) + ccon * 8.0 * pq * xfx / ((pq * pq + xfx * xfx) * (pq * pq + xfx * xfx));
            
            double gval = fpre * xrho * (2.0 * fpex1 + rho * fpexx1 * xrho - fpex1 * 7.0 / 6.0);
            
            grho_aa[i] += gval;
        }
    }
    
    void
    VWN3FuncCubicHessianAB(       CXCCubicHessianGrid& xcCubicHessianGrid,
                            const double          factor,
                            const CDensityGrid&   densityGrid)
    {
        // paramagnetic fitting factors
        
        const double spinpolf = 1.92366105093154;

        const double fourthree   = 1.333333333333333;
        
        // paramagnetic fitting factors
        
        double pa = 0.0621814, pb = 13.0720, pc = 42.7198, px0 = -0.4092860;
        
        double pq = std::sqrt(4.0 * pc - pb * pb);
        
        double pxf0 = px0 * px0 + pb * px0 + pc;
        
        double pyf0 = pq / (pb + 2.0 * px0);
        
        double b = px0 / pxf0;
        
        double c = pxf0 * pyf0;
        
        double acon= b * pb - 1.0;
        
        double bcon = 2.0 * acon + 2.0;
        
        double ccon = 2.0 * pb * (1.0 / pq - px0 / c);

        // Ferromagnetic fitting parameters

        double pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578, px0_f = -0.7432940;
        
        double pq_f = std::sqrt(4.0 * pc_f - pb_f* pb_f);
        
        double pxf0_f = px0_f * px0_f + pb_f * px0_f + pc_f;
        
        double pyf0_f = pq_f / (pb_f + 2.0 * px0_f);
        
        double b_f = px0_f / pxf0_f;
        
        double c_f = pxf0_f * pyf0_f;
        
        double acon_f = b_f * pb_f - 1.0;
        
        double bcon_f = 2.0 * acon_f + 2.0;
        
        double ccon_f = 2.0 * pb_f * (1.0 / pq_f - px0_f / c_f);

        // various prefactor
        
        double f16 = -1.0 / 6.0;
        
        double f76 = -7.0 / 6.0;
        
        double dcrs = std::pow(3.0 / (4.0 * mathconst::getPiValue()), -f16);
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();

        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_aaa = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa,xcvars::rhoa);
        
        auto grho_aab = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);
        
        auto grho_abb = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::rhob);
        
        // diamagnetic contribution
        
        #pragma omp simd aligned(rhoa, rhob, grho_aaa, grho_aab, grho_abb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double rho = rhoa[i] + rhob[i];

            double rho2 = rho*rho; 
            
            double rho3 = rho2*rho;

            double xrho  = dcrs * f16 * std::pow(rho, f76);

            double xxrho = dcrs*std::pow(rho,-13.0/6.0)*7.0/36.0;

            double zeta = (rhoa[i] - rhob[i] ) / rho;

            double spA = 2*rhob[i]/rho2; 

            double spB =-2*rhoa[i]/rho2;

            double spAA = -4*rhob[i]/rho3;

            double spAB = 2*(rho-2*rhob[i])/rho3;

            double spBB = 4*rhoa[i]/rho3;

            double zeta2  = zeta*zeta;
            
            double zeta3  = zeta2*zeta;
            
            double zeta4  = zeta3*zeta;

            // Paramagnetic 
            
            double x = dcrs * std::pow(rho, f16);

            double xf = x * x + pb * x + pc;
            
            double xfx = 2.0 * x + pb;

            double yf = pq / xfx;

            // Derivatives of epsilon

            double fpex1 = 2.0 / x + acon * xfx / xf - bcon / (x - px0) - ccon * (2.0 * yf / xfx) / (1.0 + yf * yf);

            double fpexx1 = -2.0 / (x * x) + acon * (2.0 / xf - xfx * xfx / (xf * xf))
           
                          + bcon / ((x - px0) * (x - px0)) + ccon * 8.0 * pq * xfx / ((pq * pq + xfx * xfx) * (pq * pq + xfx * xfx));

            double fpexxx1= 4.0/(x*x*x) + acon*(2.0*xfx/(xf*xf))*(xfx*xfx/xf-3.0) 
                
                        - bcon*2.0/std::pow(x - px0,3.0)
            
                        + ccon*16.0*pq*(pq*pq-3.0*xfx*xfx)/std::pow(pq*pq + xfx*xfx,3.0);

            double diff_eps_3 = 0.5*pa*xxrho*(2*fpex1 + rho*fpexx1*xrho - 7.0/6.0*fpex1)
                                
                                + 0.5*pa*xrho*(2*fpexx1*xrho 
                       
                                +fpexx1*xrho + rho*(fpexxx1*xrho*xrho+fpexx1*xxrho)
                        
                                -7.0/6.0*fpexx1*xrho);

            // Ferromagnetic

            double xf_f = x * x + pb_f * x + pc_f;
            
            double xfx_f = 2.0 * x + pb_f;

            double yf_f = pq_f / xfx_f;

            // Derivatives of epsilon

            double fpex1_f = 2.0/x + acon*xfx_f/xf_f - bcon_f/(x - px0_f) - ccon_f*(2*yf_f/xfx_f)/(1.0 + yf_f*yf_f);

            double fpexx1_f = -2.0 / (x * x) + acon_f * (2.0 / xf_f - xfx_f * xfx_f / (xf_f * xf_f))
           
                          + bcon_f / ((x - px0_f) * (x - px0_f)) + ccon_f * 8.0 * pq_f * xfx_f / ((pq_f * pq_f + xfx_f * xfx_f) 
                          
                          * (pq_f * pq_f + xfx_f * xfx_f));

            double fpexxx1_f = 4.0/(x*x*x) + acon_f*(2.0*xfx_f/(xf_f*xf_f))*(xfx_f*xfx_f/xf_f-3.0) 
                
                        - bcon_f*2.0/std::pow(x - px0_f,3.0)
            
                        + ccon_f*16.0*pq_f*(pq_f*pq_f-3.0*xfx_f*xfx_f)/std::pow(pq_f*pq_f + xfx_f*xfx_f,3.0);

            double diff_eps_3_f = 0.5*pa_f*xxrho*(2*fpex1_f + rho*fpexx1_f*xrho - 7.0/6.0*fpex1_f)
                                
                                + 0.5*pa_f*xrho*(2*fpexx1_f*xrho 
                       
                                +fpexx1_f*xrho + rho*(fpexxx1_f*xrho*xrho+fpexx1_f*xxrho)
                        
                                -7.0/6.0*fpexx1_f*xrho);

            double ef1 = fpex1_f - fpex1;

            double ef2 = fpexx1_f - fpexx1;

            double ef3 = fpexxx1_f - fpexxx1;

            double f_zeta = spinpolf * (std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf * 4.0/3.0 * (std::pow(1+zeta, 1.0/3.0)-std::pow(1-zeta, 1.0/3.0));
            
            double f_zet2 = spinpolf * 4.0/9.0 * (std::pow(1+zeta,-2.0/3.0)+std::pow(1-zeta,-2.0/3.0));
            
            double f_zet3 =-spinpolf * 8.0/27.0 * (std::pow(1+zeta,-5.0/3.0)-std::pow(1-zeta,-5.0/3.0));

            double vcfp2 = f_zeta*(ef3*zeta4+ei3*(1-zeta4));

            double eterm = f_zet1*(ef2*zeta4 + ei2*(1-zeta4)) + f_zeta*(ef2-ei2)*4*zeta3;
    
            double ef2bi = ef2-(ef1-ef0)/rho;
   
            double ei2bi = ei2-(ei1-ei0)/rho;

            double ctrm1 = 4*f_zeta*(ef2bi-ei2bi)*zeta3 +f_zet1*(ef2bi*zeta4+ei2bi*(1-zeta4));

            double ctrm2 = (ef1-ei1-ef0+ei0)*(8*f_zet1*zeta3+12*f_zeta*zeta2)
                    
                            + f_zet2*(bterm-(ef0*zeta4 + ei0*(1-zeta4)));

            double dtrm1 = f_zet2*((ef1-ef0)*zeta4+(ei1-ei0)*(1-zeta4))

                          + (8*f_zet1*zeta3+12*f_zeta*zeta2)*(ef1-ei1-ef0+ei0) +dterm/rho;
    
            double dtrm2 = ((12*f_zet2*zeta3 + 36*f_zet1*zeta2 + 24*f_zeta*zeta)*(ef0-ei0)
                        
                            + f_zet3*(ef0*zeta4+ei0*(1-zeta4)))*rho;
            
            grho_aaa[i] += (vcfp2+ eterm*spA  
                            
                            + ctrm1*(spA+spA)+ ctrm2*spA*(spA+spA) + cterm*(spAA+spAA)
                   
                            + dtrm1*(spA*spA)+ dtrm2*spA*spA*spA   + dterm*(2*spAA*spA) )*factor;
            
            grho_aab[i] += (vcfp2+ eterm*spB 
                            
                            + ctrm1*(spA+spA)+ ctrm2*spB*(spA+spA) + cterm*(spAB+spAB) 
                   
                            + dtrm1*(spA*spA)+ dtrm2*spA*spA*spB   + dterm*(2*spAB*spA) )*factor;
            
            grho_abb[i] += (vcfp2+ eterm*spA 
                            
                            + ctrm1*(spB+spB)+ ctrm2*spA*(spB+spB) + cterm*(spAB+spAB) 
                   
                            + dtrm1*(spB*spB)+ dtrm2*spB*spB*spA   + dterm*(2*spAB*spB))*factor;
        }
    }


    void
    VWN3FuncCubicHessianA(        CXCCubicHessianGrid& xcCubicHessianGrid,
                            const double          factor,
                            const CDensityGrid&   densityGrid)
    {

    }


    void
    VWN3FuncCubicHessianB(        CXCCubicHessianGrid& xcCubicHessianGrid,
                            const double          factor,
                            const CDensityGrid&   densityGrid)
    {

    }

}  // namespace vxcfuncs
