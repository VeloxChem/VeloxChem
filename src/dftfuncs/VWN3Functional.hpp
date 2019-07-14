//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef VWN3Functional_hpp
#define VWN3Functional_hpp


#include "DensityGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    /**
     Sets exchange-correlation functional to Vosko-Wilk-Nusair functional (Parameterization 3).
     
     @return the exchange-correlation functional object.
     */
    CXCFunctional setVWN3Functional();
    
    /**
     Sets primitive exchange-correlation functional to Vosko-Wilk-Nusair functional (Parameterization 3).
     
     @return the exchange-correlation functional object.
     */
    CPrimitiveFunctional setPrimitiveVWN3Functional();
    
    /**
     Implements first order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3) for dengrid::ab case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                              const double         factor,
                              const CDensityGrid&  densityGrid);
    
    
    /**
     Implements first order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3) for dengrid::lima case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncGradientA(      CXCGradientGrid& xcGradientGrid,
                             const double         factor,
                             const CDensityGrid&  densityGrid);
    
    /**
     Implements first order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3) for dengrid::lima case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncGradientB(      CXCGradientGrid& xcGradientGrid,
                           const double           factor,
                           const CDensityGrid&    densityGrid);
    
}  // namespace vxcfuncs

#endif /* VWN3Functional_hpp */
