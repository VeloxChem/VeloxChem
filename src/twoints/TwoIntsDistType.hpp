//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef TwoIntsDistType_hpp
#define TwoIntsDistType_hpp

#include <string>

/**
 Enumerate class dist2e:
 
 Defines supported two electron integrals distribution keys:
 dist2e::batch  - the batch with natural order of data
 */
enum class dist2e
{
    batch
};

/**
 Converts enumerate class value to it's string label.
 
 @param distPattern the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const dist2e distPattern)
{
    if (distPattern == dist2e::batch)
    {
        return std::string("Raw Integrals Batch");
    }
    
    return std::string("UNKNOWN");
}

#endif /* TwoIntsDistType_hpp */
