//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DummyFunctions_hpp
#define DummyFunctions_hpp

#include <vector>

#include "RecursionTerm.hpp"

namespace vlxtest {
    
    std::vector<CRecursionTerm> dummy_func_10(const CRecursionTerm& recurcionTerm);
    
    std::vector<CRecursionTerm> dummy_func_01(const CRecursionTerm& recurcionTerm);
    
    std::vector<CRecursionTerm> dummy_func_11(const CRecursionTerm& recurcionTerm);
}

#endif /* DummyFunctions_hpp */
