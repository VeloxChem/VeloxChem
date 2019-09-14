//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef XCFuncType_hpp
#define XCFuncType_hpp

#include <string>
#include "StringFormat.hpp"

/**
 Enumerate class xcfun:

 Defines supported exchange-correlation functional keys:
 xcfun::lda  - the local density approximation
 xcfun::gga  - the generalized gradient approximation
 xcfun::mgga - the meta-generalized gradient approximation
 xcfun::undefined - undefined approximation
 */
enum class xcfun
{
    lda,
    gga,
    mgga,
    undefined
};

/**
 Converts enumerate class value to its string label.

 @param functional the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const xcfun functional)
{
    if (functional == xcfun::lda) return std::string("LDA");

    if (functional == xcfun::gga) return std::string("GGA");

    if (functional == xcfun::mgga) return std::string("MGGA");

    return std::string("UNKNOWN");
}

/**
 Converts string label to its enumerate class value.

 @param label the label of enumerate class value.
 @return functional the enumerate class value.
 */
inline xcfun
to_xcfun(const std::string label)
{
    if (fstr::upcase(label) == std::string("LDA"))  return xcfun::lda;

    if (fstr::upcase(label) == std::string("GGA"))  return xcfun::gga;

    if (fstr::upcase(label) == std::string("MGGA")) return xcfun::mgga;

    return xcfun::undefined;
}

/**
 Gets number of components required for compupation of all relevant density
 contributions for specific type of exchange-correlation functions.
 
 @param xcFunctional the exchange-correlations functional type.
 @return the number of components.
 */
inline int32_t xcfun_components(const xcfun xcFunctional)
{
    if (xcFunctional == xcfun::lda) return 1;
    
    if (xcFunctional == xcfun::gga) return 4;
    
    if (xcFunctional == xcfun::mgga) return 5;
    
    return 0;
}

#endif /* XCFuncType_hpp */
