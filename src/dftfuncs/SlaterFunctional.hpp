//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef SlaterFunctional_hpp
#define SlaterFunctional_hpp

#include "DensityGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "XCFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
   
    /**
     Sets exchange-correlation functional to Slater functional.

     @return the exchange-correlation functional object.
     */
    CXCFunctional setSlaterFunctional();
    
    /**
     Sets primitive exchange-correlation functional to Slater functional.
     
     @return the exchange-correlation functional object.
     */
    CPrimitiveFunctional setPrimitiveSlaterFunctional();
    
   /**
    Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
    dengrid::ab case.

    @param xcGradientGrid the exchange-correlation gradient grid.
    @param factor the scale factor of functional contribution.
    @param densityGrid the density grid.
    */
   void SlaterFuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                             const double           factor,
                             const CDensityGrid&    densityGrid);
    
    
    /**
     Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void SlaterFuncGradientA(      CXCGradientGrid& xcGradientGrid,
                             const double           factor,
                             const CDensityGrid&    densityGrid);
    
    /**
     Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.
    
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution. 
     @param densityGrid the density grid.
    */
    void SlaterFuncGradientB(      CXCGradientGrid& xcGradientGrid,
                             const double           factor,
                             const CDensityGrid&    densityGrid);
    
    /**
    Implements second order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
    dengrid::ab case.

    @param xcHessianGrid the exchange-correlation hessian grid.
    @param factor the scale factor of functional contribution.
    @param densityGrid the density grid.
    */
   void SlaterFuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                             const double         factor,
                             const CDensityGrid&  densityGrid);
    
    
    /**
     Implements second order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void SlaterFuncHessianA(      CXCHessianGrid& xcHessianGrid,
                            const double          factor,
                            const CDensityGrid&   densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.
    
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
    */
    void SlaterFuncHessianB(      CXCHessianGrid& xcHessianGrid,
                            const double          factor,
                            const CDensityGrid&   densityGrid);
    
}  // namespace vxcfuncs

#endif /* SlaterFunctional_hpp */
