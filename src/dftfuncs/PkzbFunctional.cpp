//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "PkzbFunctional.hpp"

#include <cmath>

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
        return CPrimitiveFunctional({"PKZB"}, xcfun::mgga,
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
    PkzbFuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                         const double           factor,
                         const CDensityGrid&    densityGrid)
    {    
         // functional prefactors
            // determine number of grid points
        
        auto ngpoints = densityGrid.getNumberOfGridPoints();
        
        // set up pointers to density grid data
        
        auto rhoa = densityGrid.alphaDensity(0);
        
        auto rhob = densityGrid.betaDensity(0);

        auto grada = densityGrid.alphaDensityGradient(0);
        
        auto gradb = densityGrid.betaDensityGradient(0);
        
        auto taua = densityGrid.alphaDensityLaplacian(0);

        auto taub = densityGrid.betaDensityLaplacian(0);
        
        // set up pointers to functional data
        
        auto fexc = xcGradientGrid.xcFunctionalValues();
        
        auto grhoa = xcGradientGrid.xcGradientValues(xcvars::rhoa);
        
        auto grhob = xcGradientGrid.xcGradientValues(xcvars::rhob);
        
        auto ggrada = xcGradientGrid.xcGradientValues(xcvars::grada);
        
        auto ggradb = xcGradientGrid.xcGradientValues(xcvars::gradb);

        auto gtaua = xcGradientGrid.xcGradientValues(xcvars::lapa);
        
        auto gtaub = xcGradientGrid.xcGradientValues(xcvars::lapb);

        #pragma omp simd aligned(rhoa, rhob,grada,gradb,taua,taub, fexc, grhoa, grhob,ggrada,ggradb,gtaua,gtaub: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            // Constants 

            const double kappa = 0.804;

            const double D = 0.133;

            const double c1 = 10 / 81;

            const double c2 = 146 / 2025;

            const double c3 = -73 / 405;

            const double c4 = D + 1 / kappa * std::pow(10/81,2);

            double epsslatera = -3 / (4*mathconst::getPiValue()) *pow(3* std::pow(mathconst::getPiValue(),2),1/3) * pow(2,4/3) * std::pow(rhoa[i], 4/3);

            double epsslaterb = -3 / (4*mathconst::getPiValue()) *pow(3* std::pow(mathconst::getPiValue(),2),1/3) * pow(2,4/3) * std::pow(rhob[i], 4/3);

            double gammaa = 4 * std::pow(grhoa[i],2) / (4*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3))*std::pow(2,8/3) * std::pow(rhoa[i],8/3);

            double gammab = 4 * std::pow(grhob[i],2) / (4*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3))*std::pow(2,8/3) * std::pow(rhob[i],8/3);

            double upsilona = 3*2 * taua[i] / (2 * std::pow(3 * std::pow(mathconst::getPiValue(),2),2/3) * std::pow(2,5 / 3) * 
                            
                            std::pow(rhoa[i],5 / 3))- 9 / 20 - gammaa / 12;

            double upsilonb = 3*2 * taub[i] / (2 * std::pow(3 * std::pow(mathconst::getPiValue(),2),2/3) * std::pow(2,5 / 3) * 
                            
                            std::pow(rhob[i],5 / 3))- 9 / 20 - gammab / 12;

            double omegaa = c1 * gammaa + c2 * std::pow(upsilona,2) + c3 * upsilona * gammaa + c4 * std::pow(gammaa,2);

            double omegab = c1 * gammab + c2 * std::pow(upsilonb,2) + c3 * upsilonb * gammab + c4 * std::pow(gammab,2);

            double fa = 1 + kappa - kappa / (1+ omegaa/kappa);

            double fb = 1 + kappa - kappa / (1+ omegab/kappa);

            // derivatives

            double diffepsa = -std::pow(3*std::pow(mathconst::getPiValue(),2),1/3) / mathconst::getPiValue() * std::pow(2,4/3) * std::pow(rhoa[i],1/3);

            double diffprhoa = -8 / 3 * 4 * std::pow(grada[i],2) / (4*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3))*std::pow(2,8/3) * std::pow(rhoa[i],11/3);

            double diffqrhoa = -5 / 3 * 6 * taua[i] / (2*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3)) * std::pow(2, 5 / 3) * std::pow(rhob[i], 8 / 3) 
                            
                             - 1 / 12 * diffprhoa; 


            double diffpgrada = 8 * grada[i] / (4*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3))*std::pow(2,8/3) * std::pow(rhob[i],8/3);

            double diffqgrada = - 1 / 12 * diffpgrada;     

            double diffqtaua = 3 * 2 / (2*std::pow(3*std::pow(mathconst::getPiValue(),2),2/3))*std::pow(2,5/3) * std::pow(rhob[i],5/3);            

            double diffomegarhoa = c1 * diffprhoa + 2 * c2 * upsilona * diffqrhoa   
                                
                                 + c3 * (diffqrhoa * gammaa + upsilona * diffprhoa) + 2 * c4 * gammaa * diffprhoa;

            double diffomegagrada = c1 * diffpgrada + 2 * c2 * upsilona * diffqgrada 
                            
                                    +  c3 * (diffqgrada * gammaa + upsilona * diffpgrada) + 2 * c4 * gammaa * diffpgrada;
                     
            double diffomegataua = 2 * c2 * upsilona * diffqtaua + c3 * diffqtaua * gammaa;

            double difffomegaa = std::pow(1 + omegaa / kappa,-2);

            double difffrhoa = difffomegaa * diffomegarhoa;

            double diffgrada = difffomegaa * diffomegagrada;

            double diffftaua = difffomegaa * diffomegataua;

            fexc[i] += fa * epsslatera + fb * epsslaterb ;
            
            grhoa[i] += 0.5 * (diffepsa * fa + epsslatera * difffrhoa);
            
            grhob[i] += 0.5 * (diffepsa * fa + epsslatera * difffrhoa);

            ggrada[i] += 0.5 * epsslatera * diffgrada;

            ggradb[i] += 0.5 * epsslatera * diffgrada;

            gtaua[i] += 0.5 * epsslatera * diffftaua;

            gtaub[i] += 0.5 * epsslatera * diffftaua;    
        }

        
    }
    
    void
    PkzbFuncGradientA(      CXCGradientGrid& xcGradientGrid,
                        const double           factor,
                        const CDensityGrid&    densityGrid)
    {

    }
    
    void
    PkzbFuncGradientB(      CXCGradientGrid& xcGradientGrid,
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
    PkzbFuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {

    }
    
    void
    PkzbFuncHessianA(      CXCHessianGrid& xcHessianGrid,
                     const double          factor,
                     const CDensityGrid&   densityGrid)
    {

    }
    
    void
    PkzbFuncHessianB(      CXCHessianGrid& xcHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {
        
    }

    void
    PkzbFuncCubicHessianAB(      CXCCubicHessianGrid& xcCubicHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {

    }
    
    void
    PkzbFuncCubicHessianA(      CXCCubicHessianGrid& xcCubicHessianGrid,
                     const double          factor,
                     const CDensityGrid&   densityGrid)
    {

    }
    
    void
    PkzbFuncCubicHessianB(      CXCCubicHessianGrid& xcCubicHessianGrid,
                      const double          factor,
                      const CDensityGrid&   densityGrid)
    {
        
    }
    
}  // namespace vxcfuncs
