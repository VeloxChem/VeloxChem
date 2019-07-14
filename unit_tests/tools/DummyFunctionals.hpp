//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DummyFunctionals_hpp
#define DummyFunctionals_hpp

#include "XCGradientGrid.hpp"
#include "DensityGrid.hpp"

namespace vlxtest {
    
    void dummy_fvxc_ab(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc_a(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc_b(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
}  // namespace vlxtest


#endif /* DummyFunctionals_hpp */
