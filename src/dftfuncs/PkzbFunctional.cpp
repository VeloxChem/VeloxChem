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

#include <iostream>


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

        auto gtaua = xcGradientGrid.xcGradientValues(xcvars::taua);
        
        auto gtaub = xcGradientGrid.xcGradientValues(xcvars::taub);

        #pragma omp simd aligned(rhoa, rhob,grada,gradb,taua,taub, fexc, grhoa, grhob, ggrada, ggradb, gtaua, gtaub: VLX_ALIGN)
        for (int32_t i = 0; i < ngpoints; i++)
        {
            // Constants 
            
            const double kappa = 0.804;

            const double D = 0.133;

            const double c1 = 10.0 / 81.0;

            const double c2 = 146.0 / 2025.0;

            const double c3 = -73.0 / 405.0;

            const double c4 = D + 1.0 / kappa * std::pow(10.0/81.0,2);

            double pi_sqr = std::pow(mathconst::getPiValue(), 2.0);

            double C = -3.0 / (4.0 * mathconst::getPiValue()) *pow(3.0* pi_sqr, 1.0/3.0) * std::pow(2.0,4.0/3.0);

            double epsslatera = C *  std::pow(rhoa[i], 4.0/3.0);

            double epsslaterb = C * std::pow(rhob[i], 4.0/3.0);

            double den = (4.0*  std::pow(3.0 * pi_sqr, 2.0/3.0)) * std::pow(2,8/3);

            double gammaa = 4.0 * std::pow(grhoa[i],2.0) /den * std::pow(rhoa[i], 8.0/3.0);

            double gammab = 4.0 * std::pow(grhob[i],2.0) /den * std::pow(rhob[i],8.0/3.0);

            double den2 = (2.0 * std::pow(3.0 * pi_sqr, 2.0 /3.0 ) ) * std::pow(2.0 ,5.0 / 3.0);

            double upsilona = 3.0 * 2.0 * taua[i] / (den2 * std::pow(rhoa[i], 5.0 / 3.0))- 9.0 / 20.0 - gammaa / 12.0;

            double upsilonb = 3.0 * 2.0 * taub[i] / (den2 * std::pow(rhob[i], 5.0 / 3.0))- 9.0 / 20.0 - gammab / 12.0;

            double omegaa = c1 * gammaa + c2 * std::pow(upsilona,2.0) + c3 * upsilona * gammaa + c4 * std::pow(gammaa,2.0);

            double omegab = c1 * gammab + c2 * std::pow(upsilonb,2.0) + c3 * upsilonb * gammab + c4 * std::pow(gammab,2.0);

            double fa = 1.0 + kappa - kappa / (1.0 + omegaa/kappa);

            double fb = 1.0 + kappa - kappa / (1.0 + omegab/kappa);

            // derivatives

            double f1 = -std::pow(3.0 * pi_sqr, 1.0/3.0) / mathconst::getPiValue() * std::pow(2,4/3);

            double diffepsa =  f1 * std::pow(rhoa[i],1.0/3.0);

            double den3 =  (4*std::pow(3.0 * pi_sqr, 2.0/3.0))*std::pow(2.0,8.0/3.0);

            double diffprhoa = -8.0 / 3.0 * 4.0 * std::pow(grada[i],2.0) /( den3 * std::pow(rhoa[i],11/3) );

            double den4 =  (2.0 *std::pow(3.0 * pi_sqr, 2.0/3.0)) * std::pow(2.0, 5.0 / 3.0);

            double diffqrhoa = -5.0 / 3.0 * 6.0 * taua[i] / den4 * std::pow(rhob[i], 8.0 / 3.0) - 1.0 / 12.0 * diffprhoa;

            double diffpgrada = 8.0 * grada[i] / ( den3 * std::pow(rhob[i],8/3) ) ;

            double diffqgrada = - 1.0 / 12.0 * diffpgrada;     

            double den5 = (2*std::pow(3*std::pow(mathconst::getPiValue(), 2.0 ), 2.0/3.0))*std::pow(2.0,5.0/3.0);

            double diffqtaua = 3.0 * 2.0 /(den5  * std::pow(rhob[i],5.0/3.0) );            

            double diffomegarhoa = c1 * diffprhoa + 2.0 * c2 * upsilona * diffqrhoa   
                                
                                 + c3 * (diffqrhoa * gammaa + upsilona * diffprhoa) + 2.0 * c4 * gammaa * diffprhoa;

            double diffomegagrada = c1 * diffpgrada + 2.0 * c2 * upsilona * diffqgrada 
                            
                                    +  c3 * (diffqgrada * gammaa + upsilona * diffpgrada) + 2.0 * c4 * gammaa * diffpgrada;
                     
            double diffomegataua = 2.0 * c2 * upsilona * diffqtaua + c3 * diffqtaua * gammaa;

            double difffomegaa = std::pow(1.0 + omegaa / kappa,-2.0);

            double difffrhoa = difffomegaa * diffomegarhoa;

            double diffgrada = difffomegaa * diffomegagrada;

            double diffftaua = difffomegaa * diffomegataua;

            fexc[i] += fa * epsslatera + fb * epsslaterb;
            
            grhoa[i] += 0.5 * (diffepsa * fa + epsslatera * difffrhoa);
            
            grhob[i] += 0.5 * (diffepsa * fa + epsslatera * difffrhoa);

            ggrada[i] += 0.5 * epsslatera * diffgrada;

            ggradb[i] += 0.5 * epsslatera * diffgrada;

            gtaua[i] += 0.5 * epsslatera * diffftaua;
 
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
