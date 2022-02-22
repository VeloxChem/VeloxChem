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
    
    /**
     Implements second order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::ab case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncHessianAB(      CXCHessianGrid& xcHessianGrid,
                           const double         factor,
                           const CDensityGrid&  densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::lima case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncHessianA(      CXCHessianGrid& xcHessianGrid,
                          const double          factor,
                          const CDensityGrid&   densityGrid);
    
    /**
     Implements second order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3) for dengrid::lima case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncHessianB(      CXCHessianGrid& xcHessianGrid,
                          const double          factor,
                          const CDensityGrid&   densityGrid);


    /**
     Implements third order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::ab case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncCubicHessianAB(      CXCCubicHessianGrid& xcCubicHessianGrid,
                                const double         factor,
                                const CDensityGrid&  densityGrid);

        /**
     Implements third order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::a case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncCubicHessianA(      CXCCubicHessianGrid& xcCubicHessianGrid,
                                const double         factor,
                                const CDensityGrid&  densityGrid);
    
    /**
     Implements third order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::b case.
     
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void VWN3FuncCubicHessianB(      CXCCubicHessianGrid& xcCubicHessianGrid,
                                const double         factor,
                                const CDensityGrid&  densityGrid);
    
}  // namespace vxcfuncs

#endif /* VWN3Functional_hpp */
