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

#include "DummyFunctionals.hpp"

namespace vlxtest {
    
    void dummy_fvxc_ab(      CXCGradientGrid& xcGradientGrid,
                       const double           factor,
                       const CDensityGrid&    densityGrid)
    {
        xcGradientGrid.zero(); 
    }
    
    void dummy_fvxc_a(      CXCGradientGrid& xcGradientGrid,
                      const double           factor,
                      const CDensityGrid&    densityGrid)
    {
        xcGradientGrid.zero();
    }
    
    void dummy_fvxc_b(      CXCGradientGrid& xcGradientGrid,
                      const double           factor,
                      const CDensityGrid&    densityGrid)
    {
        xcGradientGrid.zero();
    }
    
    void dummy_fvxc2_ab(      CXCHessianGrid& xcHessianGrid,
                   const double          factor,
                   const CDensityGrid&   densityGrid)
    {
        xcHessianGrid.zero();
    }
    
    void dummy_fvxc2_a(      CXCHessianGrid& xcHessianGrid,
                       const double          factor,
                       const CDensityGrid&   densityGrid)
    {
        xcHessianGrid.zero();
    }
    
    void dummy_fvxc2_b(      CXCHessianGrid& xcHessianGrid,
                       const double           factor,
                       const CDensityGrid&    densityGrid)
    {
        xcHessianGrid.zero(); 
    }

    void dummy_fvxc3_ab(      CXCCubicHessianGrid& xcCubicHessianGrid,
                   const double          factor,
                   const CDensityGrid&   densityGrid)
    {
        xcCubicHessianGrid.zero();
    }
    
    void dummy_fvxc3_a(      CXCCubicHessianGrid& xcCubicHessianGrid,
                       const double          factor,
                       const CDensityGrid&   densityGrid)
    {
        xcCubicHessianGrid.zero();
    }
    
    void dummy_fvxc3_b(      CXCCubicHessianGrid& xcCubicHessianGrid,
                       const double           factor,
                       const CDensityGrid&    densityGrid)
    {
        xcCubicHessianGrid.zero(); 
    }
    
}  // namespace vlxtest
