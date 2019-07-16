//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "Becke88Functional.hpp"

#include <cmath>

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
                                    &vxcfuncs::Becke88FuncGradientB);
    }
    
    void
    Becke88FuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                          const double           factor,
                          const CDensityGrid&    densityGrid)
    {
        // functional prefactors
        
        double fb = 0.0042; 
        
        double fp = 4.0 / 3.0;
        
        double fre = factor * fb;
        
        // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);
        
        auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);
        
        auto ggradb = xcGradientGrid.xcGradientValues(xcvars::gradb);
        
        #pragma omp simd aligned(rhoa, rhob, grada, gradb, fexc, grhoa, grhob, ggrada, ggradb: VLX_ALIGN)
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

            fexc[i] += fre * (fea + feb);
            
            grhoa[i] += factor * grada[i] * ff1a * faa10;
            
            ggrada[i] += factor*(ffa + xa * ff1a);
            
            grhob[i] += factor * gradb[i] * ff1b * fab10;
            
            ggradb[i] += factor * (ffb + xb * ff1b);
        }
    }
    
    void
    Becke88FuncGradientA(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {
        
    }
    
    void
    Becke88FuncGradientB(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {
       
    }
    
}  // namespace vxcfuncs
