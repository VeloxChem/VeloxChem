//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef FunctionalParser_hpp
#define FunctionalParser_hpp

#include <string>

#include "XCFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace
    
    /**
     Converts exchange-correlation functional label to exchange-correlation functional object.
     
     @param xcLabel the label of exchange-correlation functional.
     @return the exchange-correlation functional object.
     */
    CXCFunctional getExchangeCorrelationFunctional(const std::string& xcLabel);
    
}  // namespace vxcfuncs

#endif /* FunctionalParser_hpp */
