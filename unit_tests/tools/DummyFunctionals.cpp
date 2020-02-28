//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
    
    void
    dummy_fvxc2_ab(      CXCHessianGrid& xcHessianGrid,
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
    
}  // namespace vlxtest
