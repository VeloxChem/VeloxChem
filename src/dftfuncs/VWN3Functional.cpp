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
        
        
        // paramagnetic fitting factors
        
        const double spinpolf = 1.92366105093154;

        const double fourthree   = 1.333333333333333;

        double px0 = -0.4092860, pa = 0.0621814, pb = 13.0720, pc = 42.7198;

        double px0_f = -0.7432940, pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578;

        // Paramagnetic 

        double Q  = std::sqrt(4.0*pc - pb*pb); 

        double XF0I = px0*px0 + pb*px0 + pc;

        double YF0I = Q/(pb + 2.0*px0); 

        double DCRS = std::pow(3.0/(4.0*mathconst::getPiValue()),1.0/6.0);

        double B    = px0/XF0I;
        
        double C    = XF0I*YF0I;

        double ACON = B*pb - 1.0; 

        double BCON = 2.0*ACON + 2.0; 

        double CCON = 2.0*pb*(1.0/Q - px0/C);

        // Ferromagnetic

        double Q_f = std::sqrt(4.0 * pc_f - pb_f*pb_f); 

        double XF0I_f = px0_f * px0_f + pb_f * px0_f + pc_f;

        double YF0I_f = Q_f/(pb_f + 2.0*px0_f); 

        double B_f    = px0_f / XF0I_f;
        
        double C_f    = XF0I_f * YF0I_f;

        double ACON_f = B_f * pb_f - 1.0; 

        double BCON_f = 2.0 * ACON_f + 2.0; 

        double CCON_f = 2.0* pb_f *(1.0 /Q_f - px0_f/C_f);
        
        
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

            double rho13 = std::pow(rho,1.0/3.0);

            double zeta   = (rhoa[i]-rhob[i])/rho;


            double f_zeta = spinpolf*(std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf*4.0/3.0 *(std::pow(1.0+zeta, 1.0/3.0)-std::pow(1.0-zeta, 1.0/3.0));


            // ep

            double x  = DCRS/std::sqrt(rho13);

            double xrho  = -DCRS*std::pow(rho,-7.0/6.0)/6.0;


            double xf   =   x * x +   pb * x + pc;

            double xf_f   = x * x + pb_f * x  + pc_f;


            double xfx  = 2.0 * x + pb;

            double xfx_f  = 2.0 * x + pb_f;


            double yf  = Q / xfx;

            double yf_f  = Q_f / xfx_f;



            double e1 = 2.0*std::log(x) + ACON*std::log(xf) - BCON*std::log(x - px0) + CCON*std::atan(yf);

            double e1_f = 2.0*std::log(x) + ACON_f*std::log(xf_f) - BCON*std::log(x - px0_f) + CCON_f*std::atan(yf_f);

            double ep_p0 = 0.5 * pa * e1;

            double ep_f0 = 0.5 * pa_f * e1_f;
            
            // ep_p1

            double ex1 = 2.0/x  + ACON*xfx/xf             - BCON/(x - px0)     - CCON*(2.0*yf/xfx)/(1.0 + yf*yf);
            
            double ex1_f = 2.0/x  + ACON_f * xfx_f / xf_f - BCON_f /(x - px0_f)- CCON_f * (2.0 * yf_f /xfx_f )/(1.0 + yf_f*yf_f);            

            double ep_p1 = 0.5*pa*(e1 + rho*ex1*xrho);
            
            double ep_f1 = 0.5*pa_f*(e1_f + rho*ex1_f*xrho);


            // Potential 

            double ef0   = ep_f0 - ep_p0;
            double ef1   = ep_f1 - ep_p1;

            double delta = f_zeta * ef0;
            double vcfp = f_zet1 * ef1;


            fexc[i] += (ep_p0 + delta)*rho*factor ;
            
            grhoa[i] +=  (ep_p1 + vcfp + delta*(1-zeta))*factor;

            grhob[i] +=  (ep_p1 + vcfp - delta*(1+zeta))*factor;

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
        
        const double spinpolf = 1.92366105093154;

        const double THREEFTHRD2 = 0.584822305543806;

        const double fourthree   = 1.333333333333333;

        double px0_i = -0.0047584, pa_i = -0.0337737, pb_i = 1.13107, pc_i = 13.0045;

        double px0 = -0.4092860, pa = 0.0621814, pb = 13.0720, pc = 42.7198;

        double px0_f = -0.7432940, pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578;

        // Paramagnetic 

        double Q  = std::sqrt(4.0*pc - pb*pb); 

        double XF0I = px0*px0 + pb*px0 + pc;

        double YF0I = Q/(pb + 2.0*px0); 

        double DCRS = std::pow(3.0/(4.0*mathconst::getPiValue()),1.0/6.0);

        double B    = px0/XF0I;
        
        double C    = XF0I*YF0I;

        double ACON = B*pb - 1.0; 

        double BCON = 2.0*ACON + 2.0; 

        double CCON = 2.0*pb*(1.0/Q - px0/C);

        // Ferromagnetic

        double Q_f = std::sqrt(4.0 * pc_f - pb_f*pb_f); 

        double XF0I_f = px0_f * px0_f + pb_f * px0_f + pc_f;

        double YF0I_f = Q_f/(pb_f + 2.0*px0_f); 

        double B_f    = px0_f / XF0I_f;
        
        double C_f    = XF0I_f * YF0I_f;

        double ACON_f = B_f * pb_f - 1.0; 

        double BCON_f = 2.0 * ACON_f + 2.0; 

        double CCON_f = 2.0* pb_f *(1.0 /Q_f - px0_f/C_f);

        // Interpolation

        double Q_i  = std::sqrt(4.0 * pc_i - pb_i * pb_i); 

        double XF0I_i = px0_i * px0_i + pb_i * px0_i + pc_i;

        double YF0I_i = Q_i /(pb_i + 2.0 * px0_i); 

        double B_i    = px0_i / XF0I_i;
        
        double C_i    = XF0I_i * YF0I_i;

        double ACON_i = B_i * pb_i - 1.0; 

        double BCON_i = 2.0 * ACON_i + 2.0; 

        double CCON_i = 2.0 * pb_i * (1.0/Q_i - px0_i / C_i);
        
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

            double rho2 = rho * rho;

            double rho3 = rho2 * rho;

            double rho4 = rho3 * rho;

            double rho13 = std::pow(rho,1.0/3.0);

            double zeta   = (rhoa[i]-rhob[i])/rho;

            double zeta2  = zeta*zeta;

            double zeta3  = zeta2*zeta;

            double zeta4  = zeta3*zeta;

            double f_zeta = spinpolf*(std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf*4.0/3.0 *(std::pow(1.0+zeta, 1.0/3.0)-std::pow(1.0-zeta, 1.0/3.0));

            double f_zet2 = spinpolf*4.0/9.0 *(std::pow(1.0+zeta,-2.0/3.0)+std::pow(1.0-zeta,-2.0/3.0));

            double f_zet3 =-spinpolf*8.0/27.0*(std::pow(1.0+zeta,-5.0/3.0)-std::pow(1.0-zeta,-5.0/3.0));

            // ep

            double x  = DCRS/std::sqrt(rho13);

            double xrho  = -DCRS*std::pow(rho,-7.0/6.0)/6.0;

            double xxrho = DCRS*std::pow(rho,-13.0/6.0)*7.0/36.0;


            double xf   =   x * x +   pb * x + pc;

            double xf_f   = x * x + pb_f * x  + pc_f;

            double xf_i   = x * x + pb_i * x + pc_i;


            double xfx  = 2.0 * x + pb;

            double xfx_f  = 2.0 * x + pb_f;

            double xfx_i  = 2.0 * x + pb_i ;



            double yf  = Q / xfx;

            double yf_f  = Q_f / xfx_f;

            double yf_i  = Q_i / xfx_i;


            double e1 = 2.0*std::log(x) + ACON*std::log(xf) - BCON*std::log(x - px0) + CCON*std::atan(yf);


            double e1_f = 2.0*std::log(x) + ACON_f*std::log(xf_f) - BCON*std::log(x - px0_f) + CCON_f*std::atan(yf_f);


            double e1_i = 2.0*std::log(x) + ACON_i*std::log(xf_i) - BCON_i*std::log(x - px0_i) + CCON_i*std::atan(yf_i);


            double ep_p0 = 0.5 * pa * e1;

            double ep_f0 = 0.5 * pa_f * e1_f;
            
            double ep_i0 = 0.5 * pa_i * e1_i;

            // ep_p1

            double ex1 = 2.0/x  + ACON*xfx/xf             - BCON/(x - px0)     - CCON*(2.0*yf/xfx)/(1.0 + yf*yf);
            
            double ex1_f = 2.0/x  + ACON_f * xfx_f / xf_f - BCON_f /(x - px0_f)- CCON_f * (2.0 * yf_f /xfx_f )/(1.0 + yf_f*yf_f);

            double ex1_i = 2.0/x  + ACON_i * xfx_i /xf_i - BCON_i /(x - px0_i)- CCON_i * (2.0*yf_i /xfx_i)/(1.0 + yf_i*yf_i);
            

            double ep_p1 = 0.5*pa*(e1 + rho*ex1*xrho);
            
            double ep_f1 = 0.5*pa_f*(e1_f + rho*ex1_f*xrho);

            double ep_i1 = 0.5 * pa_i * (e1_i + rho*ex1_i * xrho);


            // ep_p2

            double exx1 = -2.0/(x*x) + ACON*(2.0/xf - xfx*xfx/(xf*xf))
                        
                        + BCON/((x - px0)*(x - px0))
                    
                        + CCON*8.0*Q*xfx/((Q*Q + xfx*xfx)*(Q*Q + xfx*xfx));
            

            double exx1_f = -2.0/(x * x) + ACON_f * (2.0 / xf_f - xfx_f * xfx_f/(xf_f*xf_f))
                        
                        + BCON_f / ((x - px0_f) * (x - px0_f))
                    
                        + CCON_f * 8.0 * Q_f * xfx_f /((Q_f*Q_f + xfx_f*xfx_f)*(Q_f*Q_f + xfx_f*xfx_f));
            
            double exx1_i = -2.0/(x * x ) + ACON_i *(2.0/xf_i - xfx_i*xfx_i/(xf_i*xf_i))
                        
                        + BCON_i /((x - px0_i)*(x - px0_i))
                    
                        + CCON_i * 8.0 * Q_i * xfx_i /((Q_i * Q_i + xfx_i * xfx_i )*(Q_i * Q_i + xfx_i * xfx_i));
        
            
            double ep_p2 = 0.5*pa*xrho*(2.0*ex1 + rho*exx1*xrho - ex1*7.0/6.0);

            double ep_f2 = 0.5*pa_f*xrho*(2.0*ex1_f + rho*exx1_f*xrho - ex1_f*7.0/6.0);

            double ep_i2 = 0.5 * pa_i * xrho * (2.0*ex1_i + rho*exx1_i*xrho - ex1_i*7.0/6.0);


            // ep_p3

            double exxx1= 4.0/(x*x*x)
                        
                        + ACON*(2.0*xfx/(xf*xf))*(xfx*xfx/xf-3.0)
                        
                        - BCON*2.0/std::pow(x - px0,3.0)
                    
                        + CCON*16.0*Q*(Q*Q-3.0*xfx*xfx)/std::pow(Q*Q + xfx*xfx,3.0);

            double exxx1_f= 4.0 / (x * x * x)
                        
                        + ACON_f * (2.0 * xfx_f / (xf_f * xf_f))*(xfx_f * xfx_f /xf_f-3.0)
                        
                        - BCON_f *2.0 /std::pow(x - px0_f,3.0)
                    
                        + CCON_f * 16.0 * Q_f * (Q_f * Q_f-3.0 * xfx_f * xfx_f)/std::pow(Q_f * Q_f + xfx_f * xfx_f ,3.0);   

            double exxx1_i = 4.0/(x * x * x)
                        
                        + ACON_i * (2.0 * xfx_i /(xf_i * xf_i))*(xfx_i * xfx_i /xf_i-3.0)
                        
                        - BCON_i *2.0/std::pow(x - px0_i,3.0)
                    
                        + CCON_i * 16.0 * Q_i * (Q_i * Q_i -3.0*xfx_i*xfx_i )/std::pow(Q_i * Q_i + xfx_i * xfx_i,3.0);

            double ep_p3 = 0.5*pa*xxrho*(2.0*ex1 + rho*exx1*xrho - 7.0/6.0*ex1)
                        
                            + 0.5*pa*xrho*(2.0*exx1*xrho 
                       
                            +exx1*xrho + rho*(exxx1*xrho*xrho+exx1*xxrho)
                       
                            -7.0/6.0*exx1*xrho);

            double ep_f3 = 0.5 * pa_f * xxrho *(2.0*ex1_f + rho*exx1_f * xrho - 7.0/6.0 * ex1_f)
                        
                            + 0.5*pa_f*xrho*(2.0*exx1_f*xrho 
                       
                            +exx1_f * xrho + rho*(exxx1_f * xrho * xrho +exx1_f * xxrho)
                       
                            -7.0/6.0*exx1_f * xrho);
   
            double ep_i3 = 0.5 * pa_i * xxrho *(2.0 * ex1_i + rho * exx1_i * xrho - 7.0/6.0 * ex1_i)
                        
                            + 0.5 * pa_i * xrho * (2.0 * exx1_i * xrho 
                       
                            +exx1_i * xrho + rho*(exxx1_i * xrho * xrho + exx1_i * xxrho)
                       
                            -7.0/6.0 * exx1_i * xrho);

            // Potential 

            double ef0   = ep_f0 - ep_p0;
            double ef1   = ep_f1 - ep_p1;
            double ef2   = ep_f2 - ep_p2;
            double ef3   = ep_f3 - ep_p3;

            double ei0   = ep_i0*THREEFTHRD2;
            double ei1   = ep_i1*THREEFTHRD2;
            double ei2   = ep_i2*THREEFTHRD2;
            double ei3   = ep_i3*THREEFTHRD2;

            double vcfp = f_zeta*ef1;
            double delta= f_zet1*ef0;

            double vcf2 = f_zeta*ef2;
            double vap2 = f_zet1*ef1;
            double fac2 = f_zet2*ef0*rho;
            double zA   =  2*rhob[i]/rho2;
            double zB   = -2*rhoa[i]/rho2;
            double zAAr = -4*rhob[i]/rho2;
            double zABr =  2*zeta/rho;
            double zBBr =  4*rhoa[i]/rho2;

                    
            grho_aa[i] += ep_p2*factor + (vcf2 + vap2*(zA+zA) + fac2*zA*zA + delta*zAAr)*factor;
            
            grho_ab[i] += ep_p2*factor + (vcf2 + vap2*(zA+zB)  +fac2*zA*zB + delta*zABr)*factor;
            
            grho_bb[i] += ep_p2*factor + (vcf2 + vap2*(zB+zB) + fac2*zB*zB + delta*zBBr)*factor;
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

        // Paramagnetic parameters

        const double spinpolf = 1.92366105093154;

        const double THREEFTHRD2 = 0.584822305543806;

        const double fourthree   = 1.333333333333333;

        double px0_i = -0.0047584, pa_i = -0.0337737, pb_i = 1.13107, pc_i = 13.0045;

        double px0 = -0.4092860, pa = 0.0621814, pb = 13.0720, pc = 42.7198;

        double px0_f = -0.7432940, pa_f = 0.0310907, pb_f = 20.1231, pc_f = 101.578;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);

        auto df3000 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

        auto df2100 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);
        
        auto df1200 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

        #pragma omp simd aligned(rhoa, rhob,df3000,df2100,df1200: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            
            double rho = rhoa[i] + rhob[i];

            double rho2 = rho * rho;

            double rho3 = rho2 * rho;

            double rho4 = rho3 * rho;

            double rho13 = std::pow(rho,1.0/3.0);

            double zeta   = (rhoa[i]-rhob[i])/rho;

            double zeta2  = zeta*zeta;

            double zeta3  = zeta2*zeta;

            double zeta4  = zeta3*zeta;

            double f_zeta = spinpolf*(std::pow(1+zeta,fourthree)+std::pow(1-zeta,fourthree)-2.0);

            double f_zet1 = spinpolf*4.0/3.0 *(std::pow(1+zeta, 1.0/3.0)-std::pow(1-zeta, 1.0/3.0));

            double f_zet2 = spinpolf*4.0/9.0 *(std::pow(1+zeta,-2.0/3.0)+std::pow(1-zeta,-2.0/3.0));

            double f_zet3 =-spinpolf*8.0/27.0*(std::pow(1+zeta,-5.0/3.0)-std::pow(1-zeta,-5.0/3.0));

            // ep

            double Q  = std::sqrt(4.0*pc - pb*pb); 

            double XF0I = px0*px0 + pb*px0 + pc;

            double YF0I = Q/(pb + 2.0*px0); 

            double DCRS = std::pow(3.0/(4.0*mathconst::getPiValue()),1.0/6.0);

            double B    = px0/XF0I;
            
            double C    = XF0I*YF0I;

            double ACON = B*pb - 1.0; 

            double BCON = 2.0*ACON + 2.0; 

            double CCON = 2.0*pb*(1.0/Q - px0/C);

            double x  = DCRS/std::sqrt(rho13);

            double xrho  = -DCRS*std::pow(rho,-7.0/6.0)/6.0;

            double xxrho = +DCRS*std::pow(rho,-13.0/6.0)*7.0/36.0;

            double xf   = x*x + pb*x+pc;

            double xfx  = 2.0*x + pb;

            double yf  = Q/xfx;



            double e1 = 2.0*std::log(x) + ACON*std::log(xf) - BCON*std::log(x - px0) + CCON*std::atan(yf);

            double ep_p0 = 0.5 * pa * e1;

            // ep_p1

            double ex1 = 2.0/x  + ACON*xfx/xf - BCON/(x - px0)- CCON*(2.0*yf/xfx)/(1.0 + yf*yf);
            
            double ep_p1 = 0.5*pa*(e1 + rho*ex1*xrho);

            // ep_p2

            double exx1 = -2.0/(x*x) + ACON*(2.0/xf - xfx*xfx/(xf*xf))
                        
                        + BCON/((x - px0)*(x - px0))
                    
                        + CCON*8.0*Q*xfx/((Q*Q + xfx*xfx)*(Q*Q + xfx*xfx));
            
            double ep_p2 = 0.5*pa*xrho*(2.0*ex1 + rho*exx1*xrho - ex1*7.0/6.0);

            // ep_p3

            double exxx1= 4.0/(x*x*x)
                        
                        + ACON*(2.0*xfx/(xf*xf))*(xfx*xfx/xf-3.0)
                        
                        - BCON*2.0/std::pow(x - px0,3.0)
                    
                        + CCON*16.0*Q*(Q*Q-3.0*xfx*xfx)/std::pow(Q*Q + xfx*xfx,3.0);
   
            double ep_p3 = 0.5*pa*xxrho*(2.0*ex1 + rho*exx1*xrho - 7.0/6.0*ex1)
                        
                            + 0.5*pa*xrho*(2.0*exx1*xrho 
                       
                            +exx1*xrho + rho*(exxx1*xrho*xrho+exx1*xxrho)
                       
                            -7.0/6.0*exx1*xrho);

            // ep_f

            double Q_f = std::sqrt(4.0*pc_f - pb_f*pb_f); 

            double XF0I_f = px0_f * px0_f + pb_f * px0_f + pc_f;

            double YF0I_f = Q_f/(pb_f + 2.0*px0_f); 

            double DCRS_f = std::pow(3.0/(4.0*mathconst::getPiValue()),1.0/6.0);

            double x_f  = DCRS_f/std::sqrt(rho13);

            double xrho_f  = -DCRS_f*std::pow(rho,-7.0/6.0)/6.0;

            double xxrho_f = +DCRS_f*std::pow(rho,-13.0/6.0)*7.0/36.0;

            double xf_f   = x_f * x_f + pb_f * x_f  + pc_f;

            double xfx_f  = 2.0 * x_f + pb_f;

            double yf_f  = Q_f / xfx_f;

            double B_f    = px0_f / XF0I_f;
            
            double C_f    = XF0I_f * YF0I_f;

            double ACON_f = B_f * pb_f - 1.0; 

            double BCON_f = 2.0 * ACON_f + 2.0; 

            double CCON_f = 2.0* pb_f *(1.0 /Q_f - px0_f/C_f);

            double e1_f = 2.0*std::log(x_f) + ACON_f*std::log(xf_f) - BCON*std::log(x_f - px0_f) + CCON_f*std::atan(yf_f);

            double ep_f0 = 0.5 * pa_f * e1_f;

            // ep_f1

            double ex1_f = 2.0/x_f  + ACON_f * xfx_f / xf_f - BCON_f /(x_f - px0_f)- CCON_f * (2.0 * yf_f /xfx_f )/(1.0 + yf_f*yf_f);
            
            double ep_f1 = 0.5*pa_f*(e1_f + rho*ex1_f*xrho_f);

            // ep_f2

            double exx1_f = -2.0/(x_f * x_f) + ACON_f * (2.0 / xf_f - xfx_f * xfx_f/(xf_f*xf_f))
                        
                        + BCON_f / ((x_f - px0_f) * (x_f - px0_f))
                    
                        + CCON_f * 8.0 * Q_f * xfx_f /((Q_f*Q_f + xfx_f*xfx_f)*(Q_f*Q_f + xfx_f*xfx_f));
            
            double ep_f2 = 0.5*pa_f*xrho_f*(2.0*ex1_f + rho*exx1_f*xrho_f - ex1_f*7.0/6.0);

            // ep_f3

            double exxx1_f= 4.0 / (x_f * x_f * x_f)
                        
                        + ACON_f * (2.0 * xfx_f / (xf_f * xf_f))*(xfx_f * xfx_f /xf_f-3.0)
                        
                        - BCON_f *2.0 /std::pow(x_f - px0_f,3.0)
                    
                        + CCON_f * 16.0 * Q_f * (Q_f * Q_f-3.0 * xfx_f * xfx_f)/std::pow(Q_f * Q_f + xfx_f * xfx_f ,3.0);
   
            double ep_f3 = 0.5 * pa_f * xxrho_f *(2.0*ex1_f + rho*exx1_f * xrho_f - 7.0/6.0 * ex1_f)
                        
                            + 0.5*pa_f*xrho_f*(2.0*exx1_f*xrho_f 
                       
                            +exx1_f * xrho_f + rho*(exxx1_f * xrho_f * xrho_f +exx1_f * xxrho_f)
                       
                            -7.0/6.0*exx1_f * xrho_f);
            // ep_i

            double Q_i  = std::sqrt(4.0 * pc_i - pb_i * pb_i); 

            double XF0I_i = px0_i * px0_i + pb_i * px0_i + pc_i;

            double YF0I_i = Q_i /(pb_i + 2.0 * px0_i); 

            double DCRS_i = std::pow(3.0/(4.0*mathconst::getPiValue()),1.0/6.0);

            double x_i  = DCRS_i/std::sqrt(rho13);

            double xrho_i  = -DCRS_i * std::pow(rho,-7.0/6.0)/6.0;

            double xxrho_i = +DCRS_i * std::pow(rho,-13.0/6.0)*7.0/36.0;

            double xf_i   = x_i * x_i + pb_i * x_i + pc_i;

            double xfx_i  = 2.0 * x_i + pb_i ;

            double yf_i  = Q_i / xfx_i;

            double B_i    = px0_i / XF0I_i;
            
            double C_i    = XF0I_i * YF0I_i;

            double ACON_i = B_i * pb_i - 1.0; 

            double BCON_i = 2.0 * ACON_i + 2.0; 

            double CCON_i = 2.0 * pb_i * (1.0/Q_i - px0_i / C_i);

            double e1_i = 2.0*std::log(x_i) + ACON_i*std::log(xf_i) - BCON_i*std::log(x_i - px0_i) + CCON_i*std::atan(yf_i);

            double ep_i0 = 0.5 * pa_i * e1_i;

            // ep_i1

            double ex1_i = 2.0/x_i  + ACON_i * xfx_i /xf_i - BCON_i /(x_i - px0_i)- CCON_i * (2.0*yf_i /xfx_i)/(1.0 + yf_i*yf_i);
            
            double ep_i1 = 0.5 * pa_i * (e1_i + rho*ex1_i * xrho_i);

            // ep_i2

            double exx1_i = -2.0/(x_i * x_i ) + ACON_i *(2.0/xf_i - xfx_i*xfx_i/(xf_i*xf_i))
                        
                        + BCON_i /((x_i - px0_i)*(x_i - px0_i))
                    
                        + CCON_i * 8.0 * Q_i * xfx_i /((Q_i * Q_i + xfx_i * xfx_i )*(Q_i * Q_i + xfx_i * xfx_i));
            
            double ep_i2 = 0.5 * pa_i * xrho_i * (2.0*ex1_i + rho*exx1_i*xrho_i - ex1_i*7.0/6.0);

            // ep_i3

            double exxx1_i = 4.0/(x_i * x_i * x_i)
                        
                        + ACON_i * (2.0 * xfx_i /(xf_i * xf_i))*(xfx_i * xfx_i /xf_i-3.0)
                        
                        - BCON_i *2.0/std::pow(x_i - px0_i,3.0)
                    
                        + CCON_i * 16.0 * Q_i * (Q_i * Q_i -3.0*xfx_i*xfx_i )/std::pow(Q_i * Q_i + xfx_i * xfx_i,3.0);
   
            double ep_i3 = 0.5 * pa_i * xxrho_i *(2.0 * ex1_i + rho * exx1_i * xrho_i - 7.0/6.0 * ex1_i)
                        
                            + 0.5 * pa_i * xrho_i * (2.0 * exx1_i * xrho_i 
                       
                            +exx1_i * xrho_i + rho*(exxx1_i * xrho_i * xrho_i + exx1_i * xxrho_i)
                       
                            -7.0/6.0 * exx1_i * xrho_i);

            // Potential 

            double ef0   = ep_f0 - ep_p0;
            double ef1   = ep_f1 - ep_p1;
            double ef2   = ep_f2 - ep_p2;
            double ef3   = ep_f3 - ep_p3;

            double ei0   = ep_i0*THREEFTHRD2;
            double ei1   = ep_i1*THREEFTHRD2;
            double ei2   = ep_i2*THREEFTHRD2;
            double ei3   = ep_i3*THREEFTHRD2;


            double zA = 2*rhob[i]/rho2; /* =  2(1-zeta)/rho */
            
            double zB =-2*rhoa[i]/rho2; /* = -2(1+zeta)/rho */
            
            double zAA = -4*rhob[i]/rho3;
            
            double zAB = 2*(rhoa[i]-rhob[i])/rho3;
            
            double zBB = 4*rhoa[i]/rho3;
            
            double zAAA = 12*rhob[i]/rho4;
            
            double zAAB = zAAA - 4/rho3;
            
            double zBBB = -12*rhoa[i]/rho4;
            
            double zABB = zBBB + 4/rho3;

            double f0e1 = f_zeta*ef1;

            double f0e2 = f_zeta*ef2;

            double f0e3 = f_zeta*ef3;

            double f1e0 = f_zet1*ef0*rho;

            double f1e1 = f_zet1*ef1;

            double f1e2 = f_zet1*ef2;

            double f2e0 = f_zet2*ef0*rho;

            double f2e1 = f_zet2*ef1;

            double f3e0 = f_zet3*ef0*rho;


            df3000[i] += ep_p3*factor + (f0e3 + zAAA*f1e0 + (zAA + zAA + zAA)*f1e1 + (zA + zA + zA)*f1e2 +
                                (zAA*zA + zA*zAA + zAA*zA)*f2e0 + (zA*zA + zA*zA +zA*zA)*f2e1 +
                                zA*zA*zA*f3e0
                                )*factor;

            df2100[i] += ep_p3*factor + (f0e3 + zAAB*f1e0 + (zAA + zAB + zAB)*f1e1 + (zA + zA + zB)*f1e2 +
                            (zAA*zB +zAB*zA + zA*zAB)*f2e0 + (zA*zA + zA*zB + zA*zB)*f2e1 + 
                            zA*zA*zB*f3e0
                            )*factor;

            df1200[i] += ep_p3*factor +  (f0e3 + zABB*f1e0 + (zBB + zAB + zAB)*f1e1 + (zB + zB + zA)*f1e2 +
                        (zBB*zA +zAB*zB + zB*zAB)*f2e0 + (zB*zB + zB*zA + zB*zA)*f2e1 + 
                        zB*zB*zA*f3e0
                        )*factor;
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


