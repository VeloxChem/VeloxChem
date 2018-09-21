//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
    qqr
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
    
    return std::string("UNKNOWN");
}


#endif /* EriScreenerType_hpp */
