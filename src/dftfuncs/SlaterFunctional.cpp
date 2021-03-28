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

#include "SlaterFunctional.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
   
    CXCFunctional
    setSlaterFunctional()
    {
        return CXCFunctional({"Slater"}, xcfun::lda, 0.0, {setPrimitiveSlaterFunctional()}, {1.0});
    }
    
    CPrimitiveFunctional
    setPrimitiveSlaterFunctional()
    {
        return CPrimitiveFunctional({"Slater"}, xcfun::lda,
                                    &vxcfuncs::SlaterFuncGradientAB,
                                    &vxcfuncs::SlaterFuncGradientA,
                                    &vxcfuncs::SlaterFuncGradientB,
                                    &vxcfuncs::SlaterFuncHessianAB,
                                    &vxcfuncs::SlaterFuncHessianA,
                                    &vxcfuncs::SlaterFuncHessianB);
    }
    
    void
    SlaterFuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {
        // functional prefactors
        
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0);
        
        double fre = 0.75 * frg;
        
        double fp = 1.0 / 3.0;

        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        #pragma omp simd aligned(rhoa, rhob, fexc, grhoa, grhob: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double fxa = std::pow(rhoa[i], fp);
            
            double fxb = std::pow(rhob[i], fp);

            fexc[i] += fre * (rhoa[i] * fxa +  rhob[i] * fxb);
            
            grhoa[i] += frg * fxa;
            
            grhob[i] += frg * fxb;
        }
    }
    
    void
    SlaterFuncGradientA(      CXCGradientGrid& xcGradientGrid,
                        const double           factor,
                        const CDensityGrid&    densityGrid)
    {
        // functional prefactors
        
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0);
        
        double fre = 0.75 * frg;
        
        double fp = 1.0 / 3.0;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        #pragma omp simd aligned(rhob, fexc, grhob: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double fxb = std::pow(rhob[i], fp);
            
            fexc[i] += fre * rhob[i] * fxb;
            
            grhob[i] += frg * fxb;
        }
    }
    
    void
    SlaterFuncGradientB(      CXCGradientGrid& xcGradientGrid,
                        const double           factor,
                        const CDensityGrid&    densityGrid)
    {
        // functional prefactors
        
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0);
        
        double fre = 0.75 * frg;
        
        double fp = 1.0 / 3.0;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        #pragma omp simd aligned(rhoa, fexc, grhoa: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            double fxa = std::pow(rhoa[i], fp);
            
            fexc[i] += fre * rhoa[i] * fxa;
            
            grhoa[i] += frg * fxa;
        }
    }
    
    void
    SlaterFuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                        const double          factor,
                        const CDensityGrid&   densityGrid)
    {
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0) / 3.0;
        
        double fp = -2.0 / 3.0;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        #pragma omp simd aligned(rhoa, rhob, grho_aa, grho_bb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            grho_aa[i] += frg * std::pow(rhoa[i], fp);
            
            grho_bb[i] += frg * std::pow(rhob[i], fp);
        }
    }
    
    void
    SlaterFuncHessianA(      CXCHessianGrid& xcHessianGrid,
                       const double          factor,
                       const CDensityGrid&   densityGrid)
    {
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0) / 3.0;
        
        double fp = -2.0 / 3.0;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhob = densityGrid.betaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_bb = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::rhob);
        
        #pragma omp simd aligned(rhob, grho_bb: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            grho_bb[i] += frg * std::pow(rhob[i], fp);
        }
    }
    
    void
    SlaterFuncHessianB(      CXCHessianGrid& xcHessianGrid,
                       const double          factor,
                       const CDensityGrid&   densityGrid)
    {
        double frg = -factor * std::pow(6.0 / mathconst::getPiValue(), 1.0 / 3.0) / 3.0;
        
        double fp = -2.0 / 3.0;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        // set up pointers to functional data
        
        auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);
        
        #pragma omp simd aligned(rhoa, grho_aa: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            grho_aa[i] += frg * std::pow(rhoa[i], fp);
        }
    }
    
}  // namespace vxcfuncs
