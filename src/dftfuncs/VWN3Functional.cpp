//
//                           VELOXCHEM 1.0-RC
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
                                    &vxcfuncs::VWN3FuncHessianB);
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
    
}  // namespace vxcfuncs
