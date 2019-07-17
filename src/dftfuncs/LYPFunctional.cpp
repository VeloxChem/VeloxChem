//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
                                    &vxcfuncs::LYPFuncGradientB);
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
        
    }
    
    void
    LYPFuncGradientB(      CXCGradientGrid& xcGradientGrid,
                     const double           factor,
                     const CDensityGrid&    densityGrid)
    {
        
    }
    
}  // namespace vxcfuncs
