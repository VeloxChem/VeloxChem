//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef XCVarsType_hpp
#define XCVarsType_hpp

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
 xcvars::taua      - the alpha density laplacian
 xcvars::taub      - the beta density laplacian
 xcvars::undefined - the undefined gradient variable type
 */
enum class xcvars
{
    rhoa,
    rhob,
    grada,
    gradb,
    gradab,
    taua,
    taub, 
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
    
    if (xcvariable == xcvars::taua) return std::string("TAU_A");
    
    if (xcvariable == xcvars::taub) return std::string("TAU_B");
    
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
    
    if (fstr::upcase(label) == std::string("TAUA"))  return xcvars::taua;
    
    if (fstr::upcase(label) == std::string("TAUB"))  return xcvars::taub;
    
    return xcvars::undefined;
}

#endif /* XCVarsType_hpp */
