//
//                           VELOXCHEM 1.0-RC
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

#ifndef EriScreenerType_hpp
#define EriScreenerType_hpp

#include <string>

/**
 Enumerate class ericut:
 
 Defines supported two electron integrals screening keys:
 ericut::qq  - the Cauchy-Schwarz screening scheme
 ericut::qqr - the distance dependent Cauchy-Schwarz screening scheme
 */
enum class ericut
{
    qq,
    qqr,
    qqden,
    qqrden
};

/**
 Converts enumerate class value to it's string label.
 
 @param screenType the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const ericut screenType)
{
    if (screenType == ericut::qq)
    {
        return std::string("Cauchy-Schwarz (QQ)");
    }
    
    if (screenType == ericut::qqr)
    {
        return std::string("Modified Cauchy-Schwarz (QQR)");
    }
    
    if (screenType == ericut::qqden)
    {
        return std::string("Cauchy-Schwarz (QQ) + Density");
    }
    
    if (screenType == ericut::qqr)
    {
        return std::string("Modified Cauchy-Schwarz (QQR)");
    }
    
    if (screenType == ericut::qqrden)
    {
        return std::string("Modified Cauchy-Schwarz (QQR) + Density");
    }
    
    return std::string("UNKNOWN");
}


#endif /* EriScreenerType_hpp */
