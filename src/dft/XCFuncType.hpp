//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef XCFuncType_hpp
#define XCFuncType_hpp

#include <string>

/**
 Enumerate class xcfun:
 
 Defines supported exchange-correlation functional keys:
 xcfun::lda  - the local density approximation
 xcfun::gga  - the generalized gradient approximation
 xcfun::mgga - the meta-generalized gradient approximation
 */
enum class xcfun
{
    lda, gga, mgga
};

/**
 Converts enumerate class value to it's string label.

 @param functional the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const xcfun functional)
{
    if (functional == xcfun::lda) return std::string("LDA");
    
    if (functional == xcfun::gga) return std::string("GGA");
    
    if (functional == xcfun::mgga) return std::string("MGGA");
    
    return std::string("UNKNOWN"); 
}

#endif /* XCFuncType_hpp */
