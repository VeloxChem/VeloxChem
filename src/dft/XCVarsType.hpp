//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef XCVarsType_hpp
#define XCCarsType_hpp

#include <string>
#include "StringFormat.hpp"

/**
 Enumerate class xcvars:
 
 Defines supported exchange-correlation functional variable types:
 xcvars::rhoa      - the alpha density
 xcvars::rhob      - the beta density
 xcvars::grada     - the modulus of alpha density gradient
 xcvars::gradb     - the modulues of beta density gradient
 xcvars::gradab    - the product of alpha density gradient with beta density gradient
 xcvars::lapa      - the alpha density laplacian
 xcvars::lapb      - the beta density laplacian
 xcvars::undefined - the undefined gradient variable type
 */
enum class xcvars
{
    rhoa,
    rhob,
    grada,
    gradb,
    gradab,
    lapa,
    lapb, 
    undefined
};

/**
 Converts enumerate class value to its string label.
 
 @param xcvariable the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const xcvars xcvariable)
{
    if (xcvariable == xcvars::rhoa) return std::string("RHO_A");
    
    if (xcvariable == xcvars::rhob) return std::string("RHO_B");
    
    if (xcvariable == xcvars::grada) return std::string("GRAD_A");
    
    if (xcvariable == xcvars::gradb) return std::string("GRAD_B");
    
    if (xcvariable == xcvars::gradab) return std::string("GRAD_AB");
    
    if (xcvariable == xcvars::lapa) return std::string("LAP_A");
    
    if (xcvariable == xcvars::lapb) return std::string("LAP_B");
    
    return std::string("UNKNOWN");
}

/**
 Converts string label to its enumerate class value.
 
 @param label the label of enumerate class value.
 @return functional the enumerate class value.
 */
inline xcvars
to_xcvars(const std::string label)
{
    if (fstr::upcase(label) == std::string("RHOA"))  return xcvars::rhoa;
    
    if (fstr::upcase(label) == std::string("RHOB"))  return xcvars::rhob;
    
    if (fstr::upcase(label) == std::string("GRADA"))  return xcvars::grada;
    
    if (fstr::upcase(label) == std::string("GRADB"))  return xcvars::gradb;
    
    if (fstr::upcase(label) == std::string("GRADAB"))  return xcvars::gradab;
    
    if (fstr::upcase(label) == std::string("LAPA"))  return xcvars::lapa;
    
    if (fstr::upcase(label) == std::string("LAPB"))  return xcvars::lapb;
    
    return xcvars::undefined;
}

#endif /* XCVarsType_hpp */
