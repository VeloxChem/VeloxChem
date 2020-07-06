//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "FunctionalParser.hpp"

#include "StringFormat.hpp"
#include "SlaterFunctional.hpp"
#include "VWN3Functional.hpp"
#include "Becke88Functional.hpp"
#include "LYPFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    CXCFunctional
    getExchangeCorrelationFunctional(const std::string& xcLabel)
    {
        // pure spin-polarized Slater exchange functional
        
        if (fstr::upcase(xcLabel) == "SLATER") return vxcfuncs::setSlaterFunctional();
        
        // pure spin-polarized Vosko-Wilk-Nusair (Parameterization 3) correlation functional
        
        if (fstr::upcase(xcLabel) == "VWN3") return vxcfuncs::setVWN3Functional();
        
        // pure spin-polarized Becke (1988) exchange functional
        
        if (fstr::upcase(xcLabel) == "BECKE88") return vxcfuncs::setBecke88Functional();
        
        // pure spin-polarized Lee, Yang and Parr correlation functional
        
        if (fstr::upcase(xcLabel) == "LYP") return vxcfuncs::setLYPFunctional();
        
        // pure spin-polarized local density exchange-correlation functional
        
        if (fstr::upcase(xcLabel) == "SLDA")
        {
            return CXCFunctional({"SLDA"}, xcfun::lda, 0.0, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveVWN3Functional()}, {1.0, 1.0});
        }
        
        // pure spin-polarized Becke/Slater exchange functional
        
        if (fstr::upcase(xcLabel) == "B88X")
        {
            return CXCFunctional({"B88X"}, xcfun::gga, 0.0, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveBecke88Functional()}, {1.0, 1.0});
        }
        
        // pure spin-polarized BLYP exchange-correlation functional
        
        if (fstr::upcase(xcLabel) == "BLYP")
        {
            return CXCFunctional({"BLYP"}, xcfun::gga, 0.0, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional()},
                                 {1.0, 1.0, 1.0});
        }
        
        // pure spin-polarized hybrid B3LYP exchange-correlation functional
        
        if (fstr::upcase(xcLabel) == "B3LYP")
        {
            return CXCFunctional({"B3LYP"}, xcfun::gga, 0.2, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional(),
                                 setPrimitiveVWN3Functional()}, {0.8, 0.72, 0.81, 0.19});
        }
        
        // pure spin-polarized hybrid BHANDH exchange-correlation functional
        
        if (fstr::upcase(xcLabel) == "BHANDH")
        {
            return CXCFunctional({"BHANDH"}, xcfun::gga, 0.5, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveLYPFunctional()}, {0.5, 1.0});
        }
        
        // pure spin-polarized hybrid B3LYP exchange-correlation functional
        
        if (fstr::upcase(xcLabel) == "BHANDHLYP")
        {
            return CXCFunctional({"BHANDHLYP"}, xcfun::gga, 0.5, {setPrimitiveSlaterFunctional(),
                                 setPrimitiveBecke88Functional(), setPrimitiveLYPFunctional()},
                                 {0.5, 0.5, 1.0});
        }
        
        // FIX ME: add other functionals here... 
        
        return CXCFunctional();
    }
    
}  // namespace vxcfuncs
