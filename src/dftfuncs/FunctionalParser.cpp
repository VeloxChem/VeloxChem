//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "FunctionalParser.hpp"

#include "StringFormat.hpp"
#include "SlaterFunctional.hpp"
namespace vxcfuncs {  // vxcfuncs namespace
    
    CXCFunctional
    getExchangeCorrelationFunctional(const std::string& xcLabel)
    {
        // pure spin-polarized Slater exchange functional
        
        if (fstr::upcase(xcLabel) == "SLATER") return vxcfuncs::setSlaterFunctional();
        
        // FIX ME: add other functionals here... 
        
        return CXCFunctional();
    }
    
}  // namespace vxcfuncs
