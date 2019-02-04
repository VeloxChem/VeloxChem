//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MOIntsType_hpp
#define MOIntsType_hpp

#include <string>

/**
 Enumerate class moints:
 
 Defines supported types of molecular integrals:
 moints::oooo - the <oo|oo> integrals
 moints::ooov - the <oo|ov> integrals
 moints::oovv - the <oo|vv> integrals
 moints::ovov - the <ov|ov> integrals
 moints::ovvv - the <ov|vv> integrals
 moints::vvvv - the <vv|vv> integrals
 */
enum class moints
{
    oooo,
    ooov,
    oovv,
    ovov,
    ovvv,
    vvvv
};

/**
 Converts enumerate class value to it's string label.
 
 @param moIntsType the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const moints moIntsType)
{
    if (moIntsType == moints::oooo)
    {
        return std::string("<oo|oo>");
    }

    if (moIntsType == moints::ooov)
    {
        return std::string("<oo|ov>");
    }
    
    if (moIntsType == moints::oovv)
    {
        return std::string("<oo|vv>");
    }
    
    if (moIntsType == moints::ovov)
    {
        return std::string("<ov|ov>");
    }
    
    if (moIntsType == moints::ovvv)
    {
        return std::string("<ov|vv>");
    }
    
    if (moIntsType == moints::vvvv)
    {
        return std::string("<vv|vv>");
    }
    
    return std::string("UNKNOWN");
}

#endif /* MOIntsType_hpp */
