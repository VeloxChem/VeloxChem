//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef LYPFunctional_hpp
#define LYPFunctional_hpp

#include "DensityGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    /**
     Sets exchange-correlation functional to Lee, Yang and Parr functional.
     
     @return the exchange-correlation functional object.
     */
    CXCFunctional setLYPFunctional();
    
    /**
     Sets primitive exchange-correlation functional to Lee, Yang and Parr  functional.
     
     @return the exchange-correlation functional object.
     */
    CPrimitiveFunctional setPrimitiveLYPFunctional();
    
    /**
     Implements first order derivatives of Lee, Yang and Parr functional for dengrid::ab case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncGradientAB(      CXCGradientGrid& xcGradientGrid,
                           const double           factor,
                           const CDensityGrid&    densityGrid);
    
    
    /**
     Implements first order derivatives of Lee, Yang and Parr functional for dengrid::lima case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncGradientA(      CXCGradientGrid& xcGradientGrid,
                          const double           factor,
                          const CDensityGrid&    densityGrid);
    
    /**
     Implements first order derivatives of Lee, Yang and Parr functional for dengrid::lima case.
     
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncGradientB(      CXCGradientGrid& xcGradientGrid,
                          const double           factor,
                          const CDensityGrid&    densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Lee, Yang and Parr functional for dengrid::ab case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                          const double          factor,
                          const CDensityGrid&   densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Lee, Yang and Parr functional for dengrid::lima case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncHessianA(      CXCHessianGrid& xcHessianGrid,
                         const double          factor,
                         const CDensityGrid&   densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Lee, Yang and Parr functional for dengrid::lima case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void LYPFuncHessianB(      CXCHessianGrid& xcHessianGrid,
                         const double          factor,
                         const CDensityGrid&   densityGrid);
    
}  // namespace vxcfuncs

#endif /* LYPFunctional_hpp */
