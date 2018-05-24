//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef CMMType_hpp
#define CMMType_hpp

#include <string>

/**
 Enumerate class cmmtyp:
 
 Defines supported exchange-correlation functional keys:
 cmmtyp::original - the original CMM parameters model
 cmmtyp::enhanced - the enhanced CMM parameters model
 */
enum class cmmtyp
{
    original, enhanced
};

/**
 Converts enumerate class value to it's string label.

 @param cmmModel the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const cmmtyp cmmModel)
{
    if (cmmModel == cmmtyp::original)
        return std::string("Origial Parameterization");
    
    if (cmmModel == cmmtyp::enhanced)
        return std::string("Enhanced Parameterization");
    
    return std::string("Unknown Parameterization");
}

#endif /* CMMType_hpp */
