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

#ifndef PkzbFunctional_hpp
#define PkzbFunctional_hpp

#include "DensityGrid.hpp"
#include "XCFunctional.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"

namespace vxcfuncs {  // vxcfuncs namespace

    /**
     Sets exchange-correlation functional to Slater functional.

     @return the exchange-correlation functional object.
     */
    CXCFunctional setPkzbFunctional();

    /**
     Sets primitive exchange-correlation functional to Slater functional.

     @return the exchange-correlation functional object.
     */
    CPrimitiveFunctional setPrimitivePkzbFunctional();

    /**
     Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::ab case.

     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void PkzbFuncGradientAB(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);

    /**
     Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.

     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void PkzbFuncGradientA(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);

    /**
     Implements first order derivatives of spin-polarized Slater functional 1/2 [Ex(2 rho_a) + Ex(2rho_b)] with Ex(rho) = -3/4 (3/pi) (rho)^4/3 for
     dengrid::lima case.

     @param xcGradientGrid the exchange-correlation gradient grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
    */
    void PkzbFuncGradientB(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);

    void PkzbFuncHessianAB(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);

    /**
     Implements second order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3)  for dengrid::lima case.

     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void PkzbFuncHessianA(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);

    /**
     Implements second order derivatives of spin-polarized Vosko-Wilk-Nusair functional (Parameterization 3) for dengrid::lima case.

     @param xcHessianGrid the exchange-correlation hessian grid.
     @param factor the scale factor of functional contribution.
     @param densityGrid the density grid.
     */
    void PkzbFuncHessianB(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);

    /**
    Implements third order derivatives of PKZB exchange functional

    @param xcCubicHessianGrid the exchange-correlation hessian grid.
    @param factor the scale factor of functional contribution.
    @param densityGrid the density grid.
    */
    void PkzbFuncCubicHessianAB(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid);

    /**
     Implements third order derivatives of PKZB exchange functional

    @param xcCubicHessianGrid the exchange-correlation hessian grid.
    @param factor the scale factor of functional contribution.
    @param densityGrid the density grid.
    */
    void PkzbFuncCubicHessianA(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid);

    /**
    Implements third order derivatives of PKZB exchange functional

    @param xcCubicHessianGrid the exchange-correlation hessian grid.
    @param factor the scale factor of functional contribution.
    @param densityGrid the density grid.
    */
    void PkzbFuncCubicHessianB(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid);

}  // namespace vxcfuncs

#endif /* SlaterFunctional_hpp */
