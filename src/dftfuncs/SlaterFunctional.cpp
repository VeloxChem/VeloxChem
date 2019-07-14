//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
                                    &vxcfuncs::SlaterFuncGradientB);
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
    
}  // namespace vxcfuncs
