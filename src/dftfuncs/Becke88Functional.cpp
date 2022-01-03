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

#include "Becke88Functional.hpp"

#include <cmath>
#include <iostream>

#include "MathConst.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
   
    CXCFunctional
    setBecke88Functional()
    {
        return CXCFunctional({"Becke88"}, xcfun::gga, 0.0, {setPrimitiveBecke88Functional()}, {1.0});
    }
    
    CPrimitiveFunctional
    setPrimitiveBecke88Functional()
    {
        return CPrimitiveFunctional({"Becke88"}, xcfun::gga,
                                    &vxcfuncs::Becke88FuncGradientAB,
                                    &vxcfuncs::Becke88FuncGradientA,
                                    &vxcfuncs::Becke88FuncGradientB,
                                    &vxcfuncs::Becke88FuncHessianAB,
                                    &vxcfuncs::Becke88FuncHessianA,
                                    &vxcfuncs::Becke88FuncHessianB,
                                    &vxcfuncs::Becke88FuncCubicHessianAB,
                                    &vxcfuncs::Becke88FuncCubicHessianA,
                                    &vxcfuncs::Becke88FuncCubicHessianB); 
    }
    
    void
    Becke88FuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                          const double           factor,
                          const CDensityGrid&    densityGrid)
    {
        // functional prefactors
        
        double fb = 0.0042; 

        double fp = 4.0 / 3.0;
        
        double fre = -factor * fb;

        double BETA = 0.0042;

        double BETA2 = BETA*BETA;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto df1000 = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto df0100 = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        auto df0010 = xcGradientGrid.xcGradientValues(xcvars::grada);
        
        auto df0001 = xcGradientGrid.xcGradientValues(xcvars::gradb);
        
        #pragma omp simd aligned(rhoa, rhob, grada, gradb, fexc, df1000, df0100, df0010, df0001: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double ra = std::pow(rhoa[i], fp);
            
            double xa = grada[i] / ra;
            
            double faa10 = -fp * xa / rhoa[i];
            
            double rb = std::pow(rhob[i], fp);
            
            double xb = gradb[i] / rb;
            
            double fab10 = -fp * xb / rhob[i];
    
            double fda = 1.0 + 6.0 * xa * fb * std::asinh(xa);
            
            double fea = ra * xa * xa / fda;
            
            double fsqa = std::sqrt(1.0 + xa * xa);
            
            double ffa = -xa * fb / fda;
            
            double ff1a = fb * (6.0 * xa * xa * fb - fsqa) / (fsqa * fda * fda);
            
            double fdb = 1.0 + 6.0 * xb * fb * std::asinh(xb);
            
            double feb = rb * xb * xb / fdb;
            
            double fsqb = std::sqrt(1.0 + xb * xb);
            
            double ffb = -xb * fb / fdb;
            
            double ff1b = fb * (6.0 * xb * xb * fb - fsqb) / (fsqb * fdb * fdb);


            // Alpha 
            
            double alpha = grada[i]*std::pow(rhoa[i],-4.0/3.0);

            double a2 = alpha*alpha;
            
            double a3 = alpha*a2;
            
            double a4 = alpha*a3;

            double a5 = alpha*a4;

            double asha = std::asinh(alpha);

            double asha2= asha*asha;

            double alpha10 = -4.0/3.0*alpha/rhoa[i];   

            double alpha20 = -7.0/3.0*alpha10/rhoa[i];

            double alpha30 = -10./3.0*alpha20/rhoa[i]; 

            double alpha01 = std::pow(rhoa[i],-4.0/3.0);   

            double alpha11 = -4.0/3.0*alpha01/rhoa[i]; 

            double alpha21 = alpha20/grada[i];       

            double alpha10_2 = alpha10*alpha10;

            double sq1a2 = std::sqrt(1 + a2);

            double denom= 1 + 6*alpha*BETA*asha;

            double denom2 = denom*denom;

            double denom3 = denom2*denom;

            double denom4 = denom3*denom;

            double ff  = -alpha*BETA/(1.0 +6.0 *alpha*BETA*asha);
                  
            double ff1 = BETA*( 6.0 *a2*BETA - sq1a2)/(sq1a2*denom2);

            double ff2 = (6.0 *BETA2*(4.0 *alpha + a3*(3.0 - 12.0 *sq1a2*BETA) + 
                            2.0 *(std::pow(1.0 + a2 ,1.5) - 3.0 *a4*BETA)*asha))/ (std::pow(1.0 + a2,1.5)*denom3);
            
            double ff3 = (6.0 *BETA2*(6.0 + 5.0 *a2 + 2.0*a4 
	                    - 12.0 *alpha*(6.0 + 16.0 *a2 + 7.0 *a4)*BETA*asha + 
	                36.0*sq1a2*BETA*(-a2*(3.0 + 2.0*a2) + 6.0 *a5*BETA*asha - 
	                std::pow(1.0 + a2 ,2.0) * asha2) + 36.0 *a4*BETA2*(6.0 *(1.0 + a2) + (-1.0 + 2.0 *a2 )*asha2)))/(std::pow(1.0 + a2, 2.5)*denom4);


            fexc[i] += fre * (fea + feb);
            
            df1000[i] += factor*grada[i]* ff1 *alpha10;
            
            df0010[i] += factor*(ff + grada[i]*ff1*alpha01);
            
            df0100[i] += factor * gradb[i] * ff1b * fab10;
            
            df0001[i] += factor * (ffb + xb * ff1b);

            
        }
    }
    
    void
    Becke88FuncGradientA(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {
                // functional prefactors
        
        double fb = 0.0042; 

        double fp = 4.0 / 3.0;
        
        double fre = -factor * fb;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhob = densityGrid.betaDensity(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
                
        auto ggradb = xcGradientGrid.xcGradientValues(xcvars::gradb);
        
        #pragma omp simd aligned(rhob, gradb, fexc, grhob, ggradb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            
            double rb = std::pow(rhob[i], fp);
            
            double xb = gradb[i] / rb;
            
            double fab10 = -fp * xb / rhob[i];
                       
            double fdb = 1.0 + 6.0 * xb * fb * std::asinh(xb);
            
            double feb = rb * xb * xb / fdb;
            
            double fsqb = std::sqrt(1.0 + xb * xb);
            
            double ffb = -xb * fb / fdb;
            
            double ff1b = fb * (6.0 * xb * xb * fb - fsqb) / (fsqb * fdb * fdb);

            fexc[i] += fre * (feb);
            
            grhob[i] += factor * gradb[i] * ff1b * fab10;
            
            ggradb[i] += factor * (ffb + xb * ff1b);
         }
    }
    void
    Becke88FuncGradientB(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {
               // functional prefactors
        
        double fb = 0.0042; 
        
        double fp = 4.0 / 3.0;
        
        double fre = -factor * fb;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
                
        auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);
                
        #pragma omp simd aligned(rhoa, grada, fexc, grhoa, ggrada: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double ra = std::pow(rhoa[i], fp);
            
            double xa = grada[i] / ra;

            double faa10 = -fp * xa / rhoa[i];
            
            double fda = 1.0 + 6.0 * xa * fb * std::asinh(xa);
            
            double fea = ra * xa * xa / fda;
            
            double fsqa = std::sqrt(1.0 + xa * xa);
            
            double ffa = -xa * fb / fda;
            
            double ff1a = fb * (6.0 * xa * xa * fb - fsqa) / (fsqa * fda * fda);

            fexc[i] += fre * (fea);
            
            grhoa[i] += factor * grada[i] * ff1a * faa10;
            
            ggrada[i] += factor*(ffa + xa * ff1a);

        }
    }
    
    void Becke88FuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                              const double          factor,
                              const CDensityGrid&   densityGrid)
    {
        // functional prefactors
        
        double BETA = 0.0042;

        double BETA2 = BETA*BETA;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        //auto rhob = densityGrid.betaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        //auto gradb = densityGrid.betaDensityGradient(0);
        
        // set up pointers to functional data
        
        auto df2000 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        //auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        auto df1010 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);
        
        //auto gmix_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradb);
        
        auto df0020 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);
        
        //auto ggrad_bb = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradb);
        
        #pragma omp simd aligned(rhoa, grada, df2000, df1010, df0020: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            
            double alpha = grada[i]*std::pow(rhoa[i],-4.0/3.0);

            double a2 = alpha*alpha;
            
            double a3 = alpha*a2;
            
            double a4 = alpha*a3;

            double a5 = alpha*a4;

            double asha = std::asinh(alpha);

            double asha2= asha*asha;

            double alpha10 = -4.0/3.0*alpha/rhoa[i];   

            double alpha20 = -7.0/3.0*alpha10/rhoa[i];

            double alpha30 = -10./3.0*alpha20/rhoa[i]; 

            double alpha01 = std::pow(rhoa[i],-4.0/3.0);   

            double alpha11 = -4.0/3.0*alpha01/rhoa[i]; 

            double alpha21 = alpha20/grada[i];       

            double alpha10_2 = alpha10*alpha10;

            double sq1a2 = std::sqrt(1 + a2);

            double denom= 1 + 6*alpha*BETA*asha;

            double denom2 = denom*denom;

            double denom3 = denom2*denom;

            double denom4 = denom3*denom;
                  
            double ff1 = BETA*( 6.0 *a2*BETA - sq1a2)/(sq1a2*denom2);

            double ff2 = (6.0 *BETA2*(4.0 *alpha + a3*(3.0 - 12.0 *sq1a2*BETA) + 
                            2.0 *(std::pow(1.0 + a2 ,1.5) - 3.0 *a4*BETA)*asha))/ (std::pow(1.0 + a2,1.5)*denom3);
            
            double ff3 = (6.0 *BETA2*(6.0 + 5.0 *a2 + 2.0*a4 
	                    - 12.0 *alpha*(6.0 + 16.0 *a2 + 7.0 *a4)*BETA*asha + 
	                36.0*sq1a2*BETA*(-a2*(3.0 + 2.0*a2) + 6.0 *a5*BETA*asha - 
	                std::pow(1.0 + a2 ,2.0) * asha2) + 36.0 *a4*BETA2*(6.0 *(1.0 + a2) + (-1.0 + 2.0 *a2 )*asha2)))/(std::pow(1.0 + a2, 2.5)*denom4);

            df1010[i] += factor* (ff1 * alpha10 + grada[i] * (ff2 * alpha10* alpha01 + ff1 * alpha11));
            
            df2000[i] += factor * grada[i] * (ff2 * alpha10_2 + ff1 * alpha20);
            
            df0020[i] += factor * (2.0 * ff1 * alpha01 + grada[i] * (ff2 * alpha01 * alpha01));
            
            // FIX ME: Add beta part
        }
    }
    
    void Becke88FuncHessianA(      CXCHessianGrid& xcHessianGrid,
                             const double          factor,
                             const CDensityGrid&   densityGrid)
    {
        
    }
    
    void Becke88FuncHessianB(      CXCHessianGrid& xcHessianGrid,
                             const double          factor,
                             const CDensityGrid&   densityGrid)
    {
        
    }
    
    void Becke88FuncCubicHessianAB(     CXCCubicHessianGrid& xcCubicHessianGrid,
                                  const double          factor,
                                  const CDensityGrid&   densityGrid)
    {
     // functional prefactors
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
                
        auto grada = densityGrid.alphaDensityGradient(0);
                
        // set up pointers to functional data

        auto df2010 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::grada);

        auto df1020 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::grada);
        
        auto df3000 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);
        
        auto df0030 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::grada);        
                
        double BETA = 0.0042;

        double BETA2 = BETA*BETA;

        #pragma omp simd aligned(rhoa, grada, df2010, df1020,df3000,df0030: VLX_ALIGN)

        for (int32_t i = 0; i < ngpoints; i++)
        {

            // Alpha 
            
            double alpha = grada[i]*std::pow(rhoa[i],-4.0/3.0);

            double a2 = alpha*alpha;
            
            double a3 = alpha*a2;
            
            double a4 = alpha*a3;

            double a5 = alpha*a4;

            double asha = std::asinh(alpha);

            double asha2= asha*asha;

            double alpha10 = -4.0/3.0*alpha/rhoa[i];   

            double alpha20 = -7.0/3.0*alpha10/rhoa[i];

            double alpha30 = -10./3.0*alpha20/rhoa[i]; 

            double alpha01 = std::pow(rhoa[i],-4.0/3.0);   

            double alpha11 = -4.0/3.0*alpha01/rhoa[i]; 

            double alpha21 = alpha20/grada[i];       

            double alpha10_2 = alpha10*alpha10;

            double sq1a2 = std::sqrt(1 + a2);

            double denom= 1 + 6*alpha*BETA*asha;

            double denom2 = denom*denom;

            double denom3 = denom2*denom;

            double denom4 = denom3*denom;
                  
            double ff1 = BETA*( 6.0 *a2*BETA - sq1a2)/(sq1a2*denom2);

            double ff2 = (6.0 *BETA2*(4.0 *alpha + a3*(3.0 - 12.0 *sq1a2*BETA) + 
                            2.0 *(std::pow(1.0 + a2 ,1.5) - 3.0 *a4*BETA)*asha))/ (std::pow(1.0 + a2,1.5)*denom3);
            
            double ff3 = (6.0 *BETA2*(6.0 + 5.0 *a2 + 2.0*a4 
	                    - 12.0 *alpha*(6.0 + 16.0 *a2 + 7.0 *a4)*BETA*asha + 
	                36.0*sq1a2*BETA*(-a2*(3.0 + 2.0*a2) + 6.0 *a5*BETA*asha - 
	                std::pow(1.0 + a2 ,2.0) * asha2) + 36.0 *a4*BETA2*(6.0 *(1.0 + a2) + (-1.0 + 2.0 *a2 )*asha2)))/(std::pow(1.0 + a2, 2.5)*denom4);


            df2010[i] += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          grada[i]*(ff3*alpha10_2*alpha01 + 
                                 2.0*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

            df1020[i] += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          grada[i]*(ff3*alpha01*alpha01*alpha10 + 
                                 2.0*ff2*alpha11*alpha01));
    
            df3000[i] += factor*grada[i]*(ff3 * alpha10_2*alpha10 +
                                3.0*ff2*alpha10*alpha20 +
                                ff1*alpha30);

            df0030[i] += factor*(3.0*ff2*alpha01*alpha01 +
                          grada[i]*ff3*std::pow(alpha01,3.0));

        }
    }

    void Becke88FuncCubicHessianA(     CXCCubicHessianGrid& xcCubicHessianGrid,
                                  const double          factor,
                                  const CDensityGrid&   densityGrid)
    {


    }

    void Becke88FuncCubicHessianB(     CXCCubicHessianGrid& xcCubicHessianGrid,
                                  const double          factor,
                                  const CDensityGrid&   densityGrid)
    {

    }

}  // namespace vxcfuncs
