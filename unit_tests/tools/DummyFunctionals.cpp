//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DummyFunctionals.hpp"

namespace vlxtest {
    
    void dummy_fvxc_ab(CXCGradientGrid& xcGradientGrid, const CDensityGrid& densityGrid)
    {
        xcGradientGrid.zero(); 
    }
    
    void dummy_fvxc_a(CXCGradientGrid& xcGradientGrid, const CDensityGrid& densityGrid)
    {
        xcGradientGrid.zero();
    }
    
    void dummy_fvxc_b(CXCGradientGrid& xcGradientGrid, const CDensityGrid& densityGrid)
    {
        xcGradientGrid.zero();
    }
    
}  // namespace vlxtest
