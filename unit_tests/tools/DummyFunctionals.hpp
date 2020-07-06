//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef DummyFunctionals_hpp
#define DummyFunctionals_hpp

#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "DensityGrid.hpp"

namespace vlxtest {
    
    void dummy_fvxc_ab(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc_a(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc_b(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc2_ab(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc2_a(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);
    
    void dummy_fvxc2_b(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid);
    
}  // namespace vlxtest


#endif /* DummyFunctionals_hpp */
