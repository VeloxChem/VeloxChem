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
#include "VWN3Functional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    CXCFunctional
    getExchangeCorrelationFunctional(const std::string& xcLabel)
    {
        // pure spin-polarized Slater exchange functional
        
        if (fstr::upcase(xcLabel) == "SLATER") return vxcfuncs::setSlaterFunctional();
        
        // pure spin-polarized Vosko-Wilk-Nusair (Parameterization 3) functional
        
        if (fstr::upcase(xcLabel) == "VWN3") return vxcfuncs::setVWN3Functional();
        
        // pure spin-polarized local density functional
        
        if (fstr::upcase(xcLabel) == "SLDA")
        {
            return CXCFunctional({"SLDA"}, xcfun::lda, 0.0, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveVWN3Functional()}, {1.0, 1.0});
        }
        
        // FIX ME: add other functionals here... 
        
        return CXCFunctional();
    }
    
}  // namespace vxcfuncs
